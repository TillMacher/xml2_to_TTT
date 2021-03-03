import re
import pandas as pd
import csv
from ete3 import NCBITaxa
import requests
import xmltodict
import time
import PySimpleGUI as sg
from Bio import SeqIO

################################################################################
# Define functions

def get_desired_ranks(taxid, desired_ranks):
    lineage = ncbi.get_lineage(taxid)
    lineage2ranks = ncbi.get_rank(lineage)
    ranks2lineage = dict((rank, taxid) for (taxid, rank) in lineage2ranks.items())
    return {'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}

def ncbi_taxid_request(hit):
    taxonomy_list = []
    taxid = hit[1]
    results = get_desired_ranks(taxid, desired_ranks)
    taxids = [str(taxid) for taxid in list(results.values())]

    # if the taxonomy is not present
    # DO THIS
    if '<not present>' in taxids:
        for taxid in taxids:
            if taxid != '<not present>':
                url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id=" + str(taxid)
                response = requests.get(url)
                data = xmltodict.parse(response.content)
                for entry in data['eSummaryResult']['DocSum']['Item']:
                    if entry['@Name'] == 'ScientificName':
                        name = entry['#text']
                        taxonomy_list.append(name)
                time.sleep(1)
            else:
                taxonomy_list.append("")

    # if all taxonomy information is present
    # DO THIS
    else:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id=" + ','.join(taxids)
        response = requests.get(url)
        data = xmltodict.parse(response.content)
        for entry in data['eSummaryResult']['DocSum']:
            for item in entry['Item']:
                if item['@Name'] == 'ScientificName':
                    name = item['#text']
                    taxonomy_list.append(name)
    return taxonomy_list

fasta_file = "/Users/tillmacher/Desktop/MP_Projects/Projects/Sicily_fish/8_Read_tables/_data/V_cluster/V_clustering_raw_read_table_0.01_p.fasta"
results_xml = "/Users/tillmacher/Desktop/MP_Projects/Projects/Sicily_fish/8_Read_tables/_data/V_cluster/V_clustering_raw_read_table_0.01_p.xml"
results_xlsx = "/Users/tillmacher/Desktop/MP_Projects/Projects/Sicily_fish/8_Read_tables/_data/V_cluster/Sicily_fish_taxonomy_table.xlsx"

################################################################################
# Main function

def blast_taxonomy_fetcher(fasta_file, results_xml, results_xlsx):

    with open(results_xml) as myfile:
        if '<BlastXML2' not in myfile.read():
            print('Please use the xml2 format!')
        else:
            with pd.ExcelWriter(results_xlsx) as writer:
                ################################################################################
                # Convert the xml to a hit table

                f = open(results_xml)
                hit_list = []
                for line in f:
                    # collect the query name
                    if "<query-title>" in line:
                        query = re.split('>|<', line)[2]
                    # collect the query sequence length
                    if "<query-len>" in line:
                        query_len = re.split('>|<', line)[2]
                    # collect the taxonomy id
                    if "<taxid>" in line:
                        taxid = re.split('>|<', line)[2]
                    # calculate the similarity
                    if "<identity>" in line:
                        identity = re.split('>|<', line)[2]
                        diff = int(query_len) - int(identity)
                        perc_id = round(100 - diff / int(query_len) * 100, 2)
                        hit_list.append([query, taxid, perc_id, "NCBI"])

                hit_table_df = pd.DataFrame(hit_list, columns=["IDs", "Taxonomy", "Similarity", "Status"])
                hit_table_df.to_excel(writer, sheet_name='Raw hits', index=False)

                query_list = list(set(hit_table_df["IDs"].tolist()))
                sorted_hits_dict = {}

                for query in sorted(query_list):
                    prev_hit = False
                    hit_list = hit_table_df.loc[hit_table_df['IDs'] == query][["Taxonomy", "Similarity"]].values.tolist()
                    similar_hit_list = [hit_list[0]]

                    if len([list(x) for x in set(tuple(x) for x in hit_list)]) == 1:
                        sorted_hits_dict[query] = similar_hit_list

                    else:
                        for hit in hit_list:

                            if prev_hit == False:
                                prev_hit = hit
                            else:
                                if hit[1] < prev_hit[1]:
                                    sorted_hits_dict[query] = similar_hit_list
                                    break
                                else:
                                    similar_hit_list.append(hit)
                                prev_hit = hit

                similarity_filtered_list = []
                for key, values in sorted_hits_dict.items():
                    unique_values = [list(x) for x in set(tuple(x) for x in values)]
                    if len(unique_values) == 1:
                        similarity_filtered_list.append([key] + values[0])
                    else:
                        similarity_filtered_list = similarity_filtered_list + [[key] + value for value in unique_values]

                ################################################################################
                # Download the NCBI taxonomy for all hits

                ############################################################################
                ## create the progress bar window
                layout = [[sg.Text('Progress bar')],
                          [sg.ProgressBar(1000, orientation='h', size=(20, 20), key='progressbar')],
                          [sg.Cancel()]]
                window_progress_bar = sg.Window('Progress bar', layout, keep_on_top=True)
                progress_bar = window_progress_bar['progressbar']
                progress_update = 0
                progress_increase = 1000 / len(similarity_filtered_list) + 1
                ############################################################################

                ncbi = NCBITaxa()
                desired_ranks = ['phylum', 'class', 'order', 'family', 'genus', 'species']
                hit_list_2 = []

                for hit in similarity_filtered_list:
                    try:
                        taxonomy_list = ncbi_taxid_request(hit)
                    except:
                        time.sleep(3)
                        taxonomy_list = ncbi_taxid_request(hit)

                    time.sleep(1)
                    hit_list_2.append([hit[0]] + taxonomy_list + [hit[2]] + ["NCBI"])

                    ############################################################################
                    event, values = window_progress_bar.read(timeout=10)
                    if event == 'Cancel'  or event is None:
                        print('Cancel')
                        window_progress_bar.Close()
                        raise RuntimeError
                    # update bar with loop value +1 so that bar eventually reaches the maximum
                    progress_update += progress_increase
                    progress_bar.UpdateBar(progress_update)
                    ############################################################################

                window_progress_bar.Close()

                hit_table_2_df = pd.DataFrame(hit_list_2, columns=["IDs","Phylum","Class","Order","Family","Genus","Species", "Similarity", "Status"])
                hit_table_2_df.to_excel(writer, sheet_name='Taxonomy added', index=False)


                ################################################################################
                # Filter the table according to the JAMP filtering method

                hit_list_3 = []
                for hit in hit_list_2:
                    identity = hit[-2]
                    # remove empty hits first
                    if hit[1] != '':
                        if 100 >= identity >= 98:
                            if hit[1:7] != ['']*6:
                                hit_list_3.append(hit)
                        elif 98 >= identity >= 95:
                            hit[6] = ''
                            if hit[1:7] != ['']*6:
                                hit_list_3.append(hit)
                        elif 95 >= identity >= 90:
                            hit[5], hit[6] = '', ''
                            if hit[1:7] != ['']*6:
                                hit_list_3.append(hit)
                        elif 90 >= identity >= 85:
                            hit[4], hit[5], hit[6] = '', '', ''
                            if hit[1:7] != ['']*6:
                                hit_list_3.append(hit)
                        else:
                            hit[3], hit[4], hit[5], hit[6] = '', '', '', ''
                            if hit[1:7] != ['']*6:
                                hit_list_3.append(hit)

                hit_list_3 = [list(x) for x in set(tuple(x) for x in hit_list_3)]

                hit_table_3_df = pd.DataFrame(hit_list_3, columns=["IDs","Phylum","Class","Order","Family","Genus","Species", "Similarity", "Status"])
                hit_table_3_df.to_excel(writer, sheet_name='JAMP filtering', index=False)

                ################################################################################
                # Remove remaining duplicates and ambigious taxonomies

                IDs = hit_table_3_df["IDs"].values.tolist()
                duplicates = list(set([x for n, x in enumerate(IDs) if x in IDs[:n]]))
                duplicates_dict = {}

                for hit in sorted(hit_list_3):
                    if (hit[0] in duplicates and hit[0] not in duplicates_dict.keys()):
                        duplicates_dict[hit[0]] = [hit]
                    elif hit[0] in duplicates:
                        duplicates_dict[hit[0]] = duplicates_dict[hit[0]] + [hit]

                for key in duplicates_dict.keys():
                    # remove: species
                    for entry in duplicates_dict[key]:
                        entry[6] = ""
                    if len([list(x) for x in set(tuple(x) for x in duplicates_dict[key])]) == 1:
                        duplicates_dict[key] = [list(x) for x in set(tuple(x) for x in duplicates_dict[key])][0]
                    else:
                        # remove: genus
                        for entry in duplicates_dict[key]:
                            entry[5] = ""
                        if len([list(x) for x in set(tuple(x) for x in duplicates_dict[key])]) == 1:
                            duplicates_dict[key] = [list(x) for x in set(tuple(x) for x in duplicates_dict[key])][0]
                        else:
                            # remove: family
                            for entry in duplicates_dict[key]:
                                entry[4] = ""
                            if len([list(x) for x in set(tuple(x) for x in duplicates_dict[key])]) == 1:
                                duplicates_dict[key] = [list(x) for x in set(tuple(x) for x in duplicates_dict[key])][0]
                            else:
                                # remove: order
                                for entry in duplicates_dict[key]:
                                    entry[3] = ""
                                if len([list(x) for x in set(tuple(x) for x in duplicates_dict[key])]) == 1:
                                    duplicates_dict[key] = [list(x) for x in set(tuple(x) for x in duplicates_dict[key])][0]
                                else:
                                    # remove: class
                                    for entry in duplicates_dict[key]:
                                        entry[2] = ""
                                        duplicates_dict[key] = entry

                hit_list_4 = []
                for hit in hit_list_3:
                    if hit[0] in duplicates_dict.keys():
                        hit_list_4.append(duplicates_dict[hit[0]])
                    else:
                        hit_list_4.append(hit)
                hit_list_4 = [list(x) for x in set(tuple(x) for x in hit_list_4)]
                present_OTUs = [OTU[0] for OTU in hit_list_4]

                # lastly add the OTU that had no match against the NCBI database
                with open(fasta_file, "r") as handle:
                    for record in SeqIO.parse(handle, "fasta"):
                        if record.id not in present_OTUs:
                            hit_list_4.append([record.id] + ['No Match'] * 7 + ['NCBI'])

                hit_table_4_df = pd.DataFrame(hit_list_4, columns=["IDs","Phylum","Class","Order","Family","Genus","Species", "Similarity", "Status"])
                hit_table_4_df.to_excel(writer, sheet_name='TTT hit', index=False)
