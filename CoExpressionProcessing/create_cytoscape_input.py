"""
author: Geut
date: 11.10.18
creating from the tables of correlations, files for networks to read in Ccytoscape
"""
import csv
import pandas as pd
import numpy


def cyto_files_from_correlations(list_of_chaps, correlation_table_file, output_folder, is_full_network=False):
    """
    getting a file of correlation and creating files to build a cool network
    """
    corr_data = pd.DataFrame.from_csv(correlation_table_file, header=0, sep='\t', index_col=0, parse_dates=True)

    # create the file for he proteins (nodes)
    # create_node_with_metabolics_file(list_of_chaps, corr_data)

    # create the file for the interactions (bridges)
    create_network_file(list_of_chaps, corr_data, output_folder, is_full_network)


def create_network_file(chaps, corr_data, output_folder, is_full_network=False):
    """
    creates the file that will be imported as a node attributes file
    """
    fh = open(output_folder + "Network file - only positive.txt", 'w')
    file_writer = csv.writer(fh, delimiter="\t")

    # create first row:
    header_row = [consts.SOURCE_COLUMN_NAME,
                  consts.INTERACTION_COLUMN_NAME,
                  consts.TARGET_COLUMN_NAME,
                  consts.COEXPRESSION_COLUMN_NAME,
                  consts.POSITIVITY_COLUMN_NAME]
    file_writer.writerow(header_row)

    for chap in chaps:
        for protein in corr_data.columns:
            r = corr_data.at[chap, protein]
            if not is_full_network:
                if r <= 0:
                    continue

            new_network_row_data = list()
            new_network_row_data.append(chap)
            new_network_row_data.append(consts.COEXPRESSION_COLUMN_NAME)
            new_network_row_data.append(protein)
            new_network_row_data.append(abs(r))
            is_positive = r >= 0
            new_network_row_data.append(is_positive)

            # write to the file
            file_writer.writerow(new_network_row_data)

    fh.close()


def create_node_with_metabolics_file(protein_file, attributes_file, output_folder):
    """
    creates the file that will be imported as a node attributes file
    only needed if there are attributes to add. if so then convert the stream from fh to panda.
    """

    prot_data = pd.read_csv(protein_file, header=0, sep='\t', index_col=0, parse_dates=True)
    all_attrib_data = pd.read_csv(attributes_file, header=0, sep='\t', index_col=0, parse_dates=True)

    # for every prot in attrubute go to that
    for protein_name in all_attrib_data.index:
        if protein_name in prot_data.index:
            attribute_value = all_attrib_data.at[protein_name, "process"]
            if type(attribute_value) == numpy.ndarray:
                attribute_value = convert_narray_to_string(attribute_value)
            prot_data.loc[protein_name, "metabolic_process"] = attribute_value

    # save to CSV
    prot_data.to_csv(output_folder + "all_prot_attributes_less_processes.tab", sep="\t")


def convert_narray_to_string(array):
    output = ""
    for val in array:
        if not val in output:
            output += str(val) + ", "
    return output[:-2]


def create_network_for_one_protein(file_path, prot_name):

    corr_data = pd.DataFrame.from_csv(file_path, header=0, sep='\t', index_col=0, parse_dates=True)

    # filter all but this protein
    with_one_column = corr_data[[prot_name]]

    create_network_file(chaps, with_one_column, CYTO_OUTPUT_FOLDER + prot_name + ' - ', True)


class consts:
    SOURCE_COLUMN_NAME = "source"
    INTERACTION_COLUMN_NAME = "interaction"
    TARGET_COLUMN_NAME = "target"
    COEXPRESSION_COLUMN_NAME = "coexpression"
    POSITIVITY_COLUMN_NAME = "is_positive"


# ################## code for running ###################
def a_test():
    # file_path = INPUT_PERCNY_FOLDER + 'test_file.tab'
    # cancer_data = pd.DataFrame.from_csv(file_path, header=0, sep='\t', index_col=0, parse_dates=True)
    # create_network_file(["HSPD1", "HSPE1"], cancer_data, CYTO_OUTPUT_FOLDER)

    create_node_with_metabolics_file(CYTO_OUTPUT_FOLDER + "test_protein_file.txt",
                                     CYTO_OUTPUT_FOLDER + "test_processes.txt",
                                     CYTO_OUTPUT_FOLDER)


def main_func(cmnd):
    file_path = INPUT_PERCNY_FOLDER + 'Pan_CANCER_analysis.tab'
    if cmnd == "network":
        cyto_files_from_correlations(chaps, file_path, CYTO_OUTPUT_FOLDER)
    elif cmnd == "nodes":
        create_node_with_metabolics_file(CYTO_OUTPUT_FOLDER + "all_protein_file.tab",
                                         CYTO_OUTPUT_FOLDER + "protein_processes_less_categories.txt",
                                         CYTO_OUTPUT_FOLDER)
    elif cmnd == "singles":
        create_network_for_one_protein(file_path, "MTHFD2")
        create_network_for_one_protein(file_path, "SHMT2")


# ###################################
CYTO_OUTPUT_FOLDER = 'C:\\Users\\Geut\\Desktop\\'
INPUT_PERCNY_FOLDER = 'C:\\Users\\Geut\\Desktop\\Lab Analysis\\results\\Bonferroni Results no MMR1\\'
INPUT_MATCH_FOLDER = 'C:\\Users\\Geut\\Desktop\\output\\match\\'
chaps = ['CLPP', 'HSCB', 'HTRA2', 'TRAP1', 'DNAJA3', 'HSPA9', 'HSPD1', 'DNAJC19',
         'HSPE1', 'GRPEL2', 'SPG7', 'CLPX', 'YME1L1', 'LONP1', 'AFG3L2']
mini_chaps = ['HSPD1', 'HSPE1']

main_func("network")

print ' - done cytoscaping!'
