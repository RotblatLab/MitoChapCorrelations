"""
author: geut
date: 2.9.18
to make pan cancer data from the data filtered by liron's script
"""
import pandas as pd
import numpy as np

# remove the r uses. all things with *
#*from rpy2.robjects import pandas2ri
#*from rpy2.robjects.packages import importr

#*pandas2ri.activate()
#*pheatmap = importr('pheatmap')


def divide_table_to_cancer_and_normal(df):
    """

    :return:
    """
    cancer_lbls = []
    normal_lbls = []
    for index in df.index:
        if 'Cancer' in index:
            cancer_lbls.append(index)
        else:
            normal_lbls.append(index)

    c_data = df.loc[cancer_lbls]
    n_data = df.loc[normal_lbls]

    return c_data, n_data


def copy_initialized_df(df_all_tissues):
    """

    :param df_example:
    :return:
    """
    # copy the skeleton
    key = df_all_tissues.keys()[0]
    first_tissue = df_all_tissues[key]
    clean_df = first_tissue.copy()

    # init the values
    for row in clean_df.index:
        for col in clean_df.columns:
            clean_df.at[row, col] = 0

    return clean_df


def get_all_tissue_values_per_cell(dic ,row, col):
    """

    :param cancer_dic:
    :param row:
    :param col:
    :return:
    """

    value_list = []
    for tissue_name in dic.keys():
        tissue_data = dic[tissue_name]
        cell_value = tissue_data.at[row, col]
        value_list.append(cell_value)

    return value_list


def produce_pan_cancer_table(tissue_list, base_path, path_dic, output_path, analysis='percentile'):
    """
    stuff
    :return:
    """
    cancer_dic = {}
    normal_dic = {}

    # load all the cancers tables, and divide to normal and cancer dictionaries
    for tissue in tissue_list:
        path = ''
        if analysis == 'bonf':
            path = base_path + tissue + '\\' + path_dic[tissue]
        else:
            path = base_path + tissue + '\\' + path_dic[tissue]
        curr_tissue_data = pd.DataFrame.from_csv(path, header=0, sep='\t', index_col=0, parse_dates=True)
        cancer_data, normal_data = divide_table_to_cancer_and_normal(curr_tissue_data)
        cancer_dic[tissue] = cancer_data
        normal_dic[tissue] = normal_data

    # copy chap and prot to a df with them, fill with zeros
    pan_cancer_df = copy_initialized_df(cancer_dic)
    pan_normal_df = copy_initialized_df(normal_dic)

    # fill the dfs for cancer and normal
    for row in pan_cancer_df.index:
        for col in pan_cancer_df.columns:
            value_list = get_all_tissue_values_per_cell(cancer_dic, row, col)
            median_val = np.median(value_list)
            pan_cancer_df.at[row, col] = median_val

    for row in pan_normal_df.index:
        for col in pan_normal_df.columns:
            value_list = get_all_tissue_values_per_cell(normal_dic, row, col)
            median_val = np.median(value_list)
            pan_normal_df.at[row, col] = median_val

    # create tabs files
    pan_cancer_df.to_csv(output_path + 'Pan_CANCER_analysis.tab', sep='\t')
    pan_normal_df.to_csv(output_path + 'Pan_NORMAL_analysis.tab', sep='\t')


def a_test():
    path = 'C:\\Users\\Geut\\Desktop\\temp\\Results\\Thyroid_Carcinoma\\percentile\\Thyroid_Carcinoma_30_percentile_Normal58_Cancer502_corr_ATCG.tab'
    dff = pd.DataFrame.from_csv(path, header=0, sep='\t', index_col=0, parse_dates=True)
    dff = dff.iloc[:5, :7]
    print dff.to_string()
    dic = {'1': dff }
    pan_cancer_df = copy_initialized_df(dic)
    print pan_cancer_df.to_string()


def run_main():
    # run with percentile data
    all_result_path = 'C:\\Users\\Geut\\Desktop\\Lab Analysis\\Results chapXchap\\' # not up to date
    produce_pan_cancer_table(tissue_list, all_result_path, percen_paths,
                             'C:\\Users\\Geut\\Desktop\\output\\percentile\\')
    print 'first done'
    produce_pan_cancer_table(tissue_list, all_result_path, match_paths,
                             'C:\\Users\\Geut\\Desktop\\output\\match\\', 'match')
    print 'second done'


def run_single_data():
    #input_path = 'C:\\Users\\Geut\\Desktop\\Lab Analysis\\results\\Results chapXchap bonf\\'
    #produce_pan_cancer_table(tissue_list, input_path, bonffr_paths,
    #                         input_path, 'bonf')
    input_path2 = 'C:\\Users\\Geut\\Desktop\\Lab Analysis\\results\\Results Matched FDR\\'
    produce_pan_cancer_table(tissue_list, input_path2, chapchap_file_paths,
                             input_path2, 'bonf')
    print 'Panifying - done'


# region consts and strings

tissue_list = ['Breast Invasive Carcinoma',
               'Colon_Adenocarcinoma',
               'Head_and_Neck_Squamous_Cell_Carcinoma',
               'Kidney_Chromophobe',
               'Kidney_Renal_Clear_Cell_Carcinoma',
               'Kidney_Renal_Papillary_Cell_Carcinoma',
               'Liver_Hepatocellular_Carcinoma',
               'Lung_Adenocarcinoma',
               'Lung_Squamous_Cell_Carcinoma',
               'Prostate_Adenocarcinoma',
               'Stomach_Adenocarcinoma',
               'Thyroid_Carcinoma',
               'Uterine_Corpus_Endometrial_Carcinoma']

percen_paths = {'Breast Invasive Carcinoma': 'Breast Invasive Carcinoma_30_percentile_Normal113_Cancer1102_corr_ATCG.tab',
                'Colon_Adenocarcinoma': 'Colon_Adenocarcinoma_30_percentile_Normal41_Cancer478_corr_ATCG.tab',
                'Head_and_Neck_Squamous_Cell_Carcinoma': 'Head_and_Neck_Squamous_Cell_Carcinoma_30_percentile_Normal44_Cancer500_corr_ATCG.tab',
                'Kidney_Chromophobe': 'Kidney_Chromophobe_30_percentile_Normal24_Cancer65_corr_ATCG.tab',
                'Kidney_Renal_Clear_Cell_Carcinoma': 'Kidney_Renal_Clear_Cell_Carcinoma_30_percentile_Normal72_Cancer538_corr_ATCG.tab',
                'Kidney_Renal_Papillary_Cell_Carcinoma': 'Kidney_Renal_Papillary_Cell_Carcinoma_30_percentile_Normal32_Cancer288_corr_ATCG.tab',
                'Liver_Hepatocellular_Carcinoma': 'Liver_Hepatocellular_Carcinoma_30_percentile_Normal50_Cancer371_corr_ATCG.tab',
                'Lung_Adenocarcinoma': 'Lung_Adenocarcinoma_30_percentile_Normal59_Cancer533_corr_ATCG.tab',
                'Lung_Squamous_Cell_Carcinoma': 'Lung_Squamous_Cell_Carcinoma_30_percentile_Normal49_Cancer502_corr_ATCG.tab',
                'Prostate_Adenocarcinoma': 'Prostate_Adenocarcinoma_30_percentile_Normal52_Cancer498_corr_ATCG.tab',
                'Stomach_Adenocarcinoma': 'Stomach_Adenocarcinoma_30_percentile_Normal32_Cancer375_corr_ATCG.tab',
                'Thyroid_Carcinoma': 'Thyroid_Carcinoma_30_percentile_Normal58_Cancer502_corr_ATCG.tab',
                'Uterine_Corpus_Endometrial_Carcinoma': 'Uterine_Corpus_Endometrial_Carcinoma_30_percentile_Normal35_Cancer551_corr_ATCG.tab'}

bonffr_paths = {'Breast Invasive Carcinoma': 'Breast Invasive Carcinoma_-1_percentile_Normal113_Cancer1102_corr_ATCG.tab',
                'Colon_Adenocarcinoma': 'Colon_Adenocarcinoma_-1_percentile_Normal41_Cancer478_corr_ATCG.tab',
                'Head_and_Neck_Squamous_Cell_Carcinoma': 'Head_and_Neck_Squamous_Cell_Carcinoma_-1_percentile_Normal44_Cancer500_corr_ATCG.tab',
                'Kidney_Chromophobe': 'Kidney_Chromophobe_-1_percentile_Normal24_Cancer65_corr_ATCG.tab',
                'Kidney_Renal_Clear_Cell_Carcinoma': 'Kidney_Renal_Clear_Cell_Carcinoma_-1_percentile_Normal72_Cancer538_corr_ATCG.tab',
                'Kidney_Renal_Papillary_Cell_Carcinoma': 'Kidney_Renal_Papillary_Cell_Carcinoma_-1_percentile_Normal32_Cancer288_corr_ATCG.tab',
                'Liver_Hepatocellular_Carcinoma': 'Liver_Hepatocellular_Carcinoma_-1_percentile_Normal50_Cancer371_corr_ATCG.tab',
                'Lung_Adenocarcinoma': 'Lung_Adenocarcinoma_-1_percentile_Normal59_Cancer533_corr_ATCG.tab',
                'Lung_Squamous_Cell_Carcinoma': 'Lung_Squamous_Cell_Carcinoma_-1_percentile_Normal49_Cancer502_corr_ATCG.tab',
                'Prostate_Adenocarcinoma': 'Prostate_Adenocarcinoma_-1_percentile_Normal52_Cancer498_corr_ATCG.tab',
                'Stomach_Adenocarcinoma': 'Stomach_Adenocarcinoma_-1_percentile_Normal32_Cancer375_corr_ATCG.tab',
                'Thyroid_Carcinoma': 'Thyroid_Carcinoma_-1_percentile_Normal58_Cancer502_corr_ATCG.tab',
                'Uterine_Corpus_Endometrial_Carcinoma': 'Uterine_Corpus_Endometrial_Carcinoma_-1_percentile_Normal35_Cancer551_corr_ATCG.tab'}

matched_paths = {'Breast Invasive Carcinoma':             'Breast Invasive Carcinoma_matched_corr_ATCG.tab',
                'Colon_Adenocarcinoma':                  'Colon_Adenocarcinoma_matched_corr_ATCG.tab',
                'Head_and_Neck_Squamous_Cell_Carcinoma': 'Head_and_Neck_Squamous_Cell_Carcinoma_matched_corr_ATCG.tab',
                'Kidney_Chromophobe':                    'Kidney_Chromophobe_matched_corr_ATCG.tab',
                'Kidney_Renal_Clear_Cell_Carcinoma':     'Kidney_Renal_Clear_Cell_Carcinoma_matched_corr_ATCG.tab',
                'Kidney_Renal_Papillary_Cell_Carcinoma': 'Kidney_Renal_Papillary_Cell_Carcinoma_matched_corr_ATCG.tab',
                'Liver_Hepatocellular_Carcinoma':        'Liver_Hepatocellular_Carcinoma_matched_corr_ATCG.tab',
                'Lung_Adenocarcinoma':                   'Lung_Adenocarcinoma_matched_corr_ATCG.tab',
                'Lung_Squamous_Cell_Carcinoma':          'Lung_Squamous_Cell_Carcinoma_matched_corr_ATCG.tab',
                'Prostate_Adenocarcinoma':               'Prostate_Adenocarcinoma_matched_corr_ATCG.tab',
                'Stomach_Adenocarcinoma':                'Stomach_Adenocarcinoma_matched_corr_ATCG.tab',
                'Thyroid_Carcinoma':                     'Thyroid_Carcinoma_matched_corr_ATCG.tab',
                'Uterine_Corpus_Endometrial_Carcinoma':  'Uterine_Corpus_Endometrial_Carcinoma_matched_corr_ATCG.tab'}

chapchap_file_paths = {}
for tissue in tissue_list:
    chapchap_file_paths[tissue] = tissue + '_match_corr_ATCG.tab'

# endregion

#run_main()
run_single_data()
print ' - done!'
