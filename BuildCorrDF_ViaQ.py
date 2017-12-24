###################################################################################################
# File: BuildCorrDF_ViaQ.py
# Version: 0.0
# Date: 30.8.2017
# Noam Bar, noam.bar@weizmann.ac.il
#
# 
# Python version: 2.7
###################################################################################################
import pickle
from queue.qp import qp
from datetime import datetime
from scipy.stats.stats import spearmanr, pearsonr
import pandas as pd
from pandas import concat
import os
import numpy as np
import argparse
from addloglevels import sethandlers
import logging
import Utils
from scipy.spatial.distance import braycurtis

log_ = logging.getLogger(__name__)
del logging

METABOLON_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Metabolon/'
LEGAL_METHODS = {'spearman':spearmanr, 'pearson':pearsonr, 'braycurtis':braycurtis}

def calc_pairwise_corr(A, B, genePos, output_dir, minimal_corr_pair, directional_pval, 
                       write_to_csv = False, write_corr = False, read_B_from_file = None):
    """
    Calculate pairwise correlation between two dataframes
    :param A: A pandas DataFrame of the form (SID, metabolites)
    :param B: A pandas DataFrame of the form (SID, EM genes)
    :param genePos: An integer indicating the starting index of the EM Genes matrix 
    :return: A dictionary with two keys:
    'corr' - correlation dataframe
    'pval' - p-value dataframe
    """
    if read_B_from_file is not None:
        print (now() + " - calc_pairwise_corr - reading B from file")
        B = Utils.Load(read_B_from_file)
    
    assert A.shape[0] == B.shape[0], "Number of rows is different."
    print (now() + " - calc_pairwise_corr - " + str(genePos))
    print (now() + " - A.shape = "  + str(A.shape))
    print (now() + " - B.shape = "  + str(B.shape))
    # keep only genes with more than one unique value (not all were null)
    valid_bac_genes = [i for i in range(B.shape[1]) if len(B.iloc[:,i].unique()) > 1]
    B = B.iloc[:, valid_bac_genes]
    print (now() + " - B.shape = "  + str(B.shape))
    resCorr = pd.DataFrame(0, index=A.columns, columns=B.columns, dtype=np.float64)
    resPval = pd.DataFrame(0, index=A.columns, columns=B.columns, dtype=np.float64)
    
    C = concat([A, B], axis=1)
    A = C.iloc[:,0:A.shape[1]].values
    B = C.iloc[:,A.shape[1]:].values

    for i in range(A.shape[1]):
        a = A[:,i]
        log_.info(str(i))
        corrs = []
        pvals = []
        for j in range(B.shape[1]):
            b = B[:,j]
            corr = None
            pval = None
            if j % 10000 == 0:
                log_.info(str(i) + ' ' + str(j))
            c = ~np.isnan(a) & ~np.isnan(b)
            if c.sum() >= minimal_corr_pair:
                corr, pval = spearmanr(np.stack((a[c],b[c]),axis=1))
            corrs.append(corr)
            pvals.append(pval)
        resCorr.iloc[i,:] = corrs
        resPval.iloc[i,:] = pvals   
        
    if directional_pval:
        resPval = resPval.apply(np.log10)
        resPval = resPval * -resCorr
    resCorr = resCorr.dropna(how='all', axis=1)
    resPval = resPval.dropna(how='all', axis=1)
    print (now() + " - resPval.shape = "  + str(resPval.shape))
    res_dict = {'corr':resCorr, 'pval':resPval}
    if write_to_csv:
        log_.info(" - writing df " + str(genePos) + " to csv")
        if write_corr:
            resCorr.T.to_csv(output_dir + 'corrs_' + str(genePos) + '.csv')
        resPval.T.to_csv(output_dir + 'pvals_' + str(genePos) + '.csv')
    else:
        log_.info(" - writing df " + str(genePos) + " to pickle")
        with open(output_dir + 'corr_pval_dict_' + str(genePos) + '.df', 'wb') as handle:
            pickle.dump(res_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    log_.info(" - one more job finished")
    return res_dict

def upload_jobs(q, metabs_df, EMGenes_df, output_dir, minimal_corr_pair, directional_pval, 
                n_cols_per_job, write_corr):
    """
    Upload jobs to the qp, breaks the EM genes matrix into sub matrices
    :param q: a qp object
    :param metabs_df: A pandas DataFrame of the form (SID, metabolites)
    :param EMGenes_df: A pandas DataFrame of the form (SID, EM genes)
    :return: A list of dictionaries with two keys:
    'corr' - correlation dataframe
    'pval' - p-value dataframe
    """
    waiton = []
    for genePos in range(0, EMGenes_df.shape[1], n_cols_per_job):
        print (now() + " - Append job #" + str(genePos))
        waiton.append(q.method(calc_pairwise_corr, 
                               (metabs_df, 
                                EMGenes_df.iloc[:,genePos:genePos+n_cols_per_job], 
                                genePos,
                                output_dir,
                                minimal_corr_pair,
                                directional_pval,
                                write_corr)))

    print (now() + " - Waiting for results")
    res = q.waitforresults(waiton)
    print (now() + " - Results are back")
    return res

def find_free_idx(idx, num_of_shuffles):
    if len(idx) == 0:
        return range(num_of_shuffles)
    free_idx = []
    shuffles_left = num_of_shuffles
    for i in range(idx[-1]):
        if i not in idx:
            free_idx.append(i)
            shuffles_left -= 1
        if shuffles_left == 0:
            return free_idx
    free_idx += range(idx[-1] + 1, idx[-1] + shuffles_left + 1)
    return free_idx

def upload_jobs_shuffled(q, metabs_df, df, output_dir, minimal_corr_pair, directional_pval, 
                         num_of_shuffels, write_corr, pass_path_to_q):
    """
    Upload jobs to the qp, breaks the EM genes matrix into sub matrices
    :param q: a qp object
    :param metabs_df: A pandas DataFrame of the form (SID, metabolites)
    :param EMGenes_df: A pandas DataFrame of the form (SID, EM genes)
    :return: A list of dictionaries with two keys:
    'corr' - correlation dataframe
    'pval' - p-value dataframe
    """
    waiton = []
#     start_idx = 0
    pval_files = os.listdir(output_dir)
    pval_files = [int((p.split('pvals_')[1]).split('.csv')[0]) for p in pval_files if p[0] == 'p']
    pval_files.sort()
    free_idx = find_free_idx(pval_files,num_of_shuffels)
#     if len(pval_files) > 0:
#         start_idx = pval_files[-1] + 1
        
    big_df_path = None
    if pass_path_to_q:
        big_df_path = output_dir + '/temp_big_df.df'
        Utils.Write(big_df_path, df)
        df = None
    
    for i in free_idx:
        print (now() + " - Append job #" + str(i))
        metabs_temp = metabs_df.copy()
        # shuffle the rows indexes
        metabs_temp.index = np.random.permutation(metabs_df.index)
        waiton.append(q.method(calc_pairwise_corr, 
                               (metabs_temp, 
                                df,
                                i,
                                output_dir,
                                minimal_corr_pair,
                                directional_pval,
                                True,
                                write_corr,
                                big_df_path)))

    print (now() + " - Waiting for results")
    res = q.wait(waiton)
    if pass_path_to_q:
        os.remove(big_df_path)
    print (now() + " - Results are back")
    return res

def concat_temp_corr_dfs(output_dir, n_cols_per_job, write_corr = False):
    """
    Read correlation/p-value dataframe and writes it to one csv file
    """
    if not os.path.exists(output_dir):
        raise output_dir + " Does not exist"
    print (now() + " - concat_temp_corr_dfs")
    if os.path.exists(output_dir + 'corr_pval_dict_0.df'):
        # read first file
        with open(output_dir + 'corr_pval_dict_0.df', 'rb') as handle:
            d = pickle.load(handle)
            if write_corr:
                d['corr'].T.to_csv(output_dir + 'corrs.csv')
            d['pval'].T.to_csv(output_dir + 'pvals.csv')
            os.remove(output_dir + 'corr_pval_dict_0.df')

    temp_files = os.listdir(output_dir)
    temp_files = [f for f in temp_files if f.split('_')[0] == 'corr']
    if len(temp_files) == 0:
        return
    files_idx = [int(f.split('corr_pval_dict_')[1].split('.df')[0]) for f in temp_files]
    files_idx.sort()
#     files_idx = files_idx[1:]
    check_list = range(n_cols_per_job, files_idx[-1]+1, n_cols_per_job)
    assert files_idx == check_list, "It appears that a temp file/s is missing"
    print (now() + " - reading and writing temp files...")
    for i in files_idx:
        if (i % 50000 == 0):
            print (now() + " - index: " + str(i))
        with open(output_dir + 'corr_pval_dict_' + str(i) + '.df', 'rb') as handle:
            d = pickle.load(handle)
        if write_corr:
            with open(output_dir + 'corrs.csv', 'a') as handle:
                d['corr'].T.to_csv(handle, header=False)
        with open(output_dir + 'pvals.csv', 'a') as handle:
            d['pval'].T.to_csv(handle, header=False)
        os.remove(output_dir + 'corr_pval_dict_' + str(i) + '.df')
    return

def remove_non_kegg_metabolites(ppmet, metabsPath):
    with open(metabsPath, 'rb') as handle:
        metabs = pickle.load(handle)
    metabs_w_kegg = [i for i in ppmet.columns if str(metabs.loc[i].KEGG) != 'nan']
    return ppmet[metabs_w_kegg]

def remove_highly_correlated_metabolites(ppmet):
    # check correlations
    # for each cluster of highly correlated metabolites keep only one
    # use Tal's code
    return ppmet

def filter_metabolites_from_file(ppmet, met_file):
    with open(met_file, 'r') as handle:
        metabs = handle.readlines()
    metabs = [int(float(c.split('\n')[0])) for c in metabs]
    return ppmet[metabs]

def now():
    return str(datetime.now())
        
def main():
    """
    Compute all correlations between EM genes matrix and metabolites and write to disk.
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-output_dir', help='Path to output directory', type=str, default=None)
    parser.add_argument('analysisType', help='What type of analysis to perform: single or shuffle', type=str, default='single')
    parser.add_argument('-data', help='What type of data: EMGenes or KO', type=str, default='EMGenes')
    parser.add_argument('-metabsPath', help='Path to metabolites file', type=str, default=METABOLON_DIR + '/tmp_files/metabs.df')
    parser.add_argument('-ppmetPath', help='Path to metabolite values file', type=str, default=METABOLON_DIR + '/tmp_files/ppmet_483x820.df')
    parser.add_argument('-EMGenesPath', help='Path to EMGenes dataframe', type=str, default=None)
    parser.add_argument('-ko', '--ko_path', help='Path to KOxsample data frame', type=str, default=METABOLON_DIR + 'KO_EMGenes.df')
    parser.add_argument('-f', '--from_file', help='Path to external file holding metabolites IDs', type=str, default=None)
    parser.add_argument('-rnk', '--remove_nonKEGG_metabs', type=bool, help='Whether to remove metabolits without KEGG IDs', default=True)
    parser.add_argument('-rhc', '--remove_high_corr_metabs', type=bool, help='Whether to remove highly correlated metabolits', default=False)
    parser.add_argument('-nc', '--n_cols_per_job', help='Number of columns per job', type=int, default=10000)
    parser.add_argument('-cp', '--minimal_corr_pair', help='Minimal number of pairs for test', type=int, default=10)
    parser.add_argument('-dp', '--directional_pval', type=bool, help='Whether to compute directional p-value', default=True)
    parser.add_argument('-oc', '--only_concat', type=bool, help='Whether to only run the concat function', default=False)
    parser.add_argument('-ns', '--num_of_shuffles', help='Number of shuffles to perform', type=int, default=100)
    parser.add_argument('-wc', '--write_corr', help='Whether to write the correlation files', type=bool, default=False)
    parser.add_argument('-pq', '--pass_path_to_q', help='Whether to pass the big df path to q', type=bool, default=True)
    parser.add_argument('-mem', '--memory_per_job', help='Amount of memory to assign each job', type=str, default='10G')
    command_args = parser.parse_args()

    if command_args.analysisType != 'single' and command_args.analysisType != 'shuffle':
        print(now() + " - [main] analysisType must be either single or shuffle.")
        return
    if command_args.data != 'EMGenes' and command_args.data != 'KO':
        print(now() + " - [main] data must be either EMGenes or KO.")
        return
    if not (os.path.exists(command_args.metabsPath) or os.path.exists(command_args.ppmetPath)):
        print(now() + " - [main] metabsPath or ppmetPath doesn't exist.")
        return
    
    if command_args.only_concat:
        concat_temp_corr_dfs(command_args.output_dir, command_args.n_cols_per_job, 
                             command_args.write_corr)
        return
    
#     return #TODO: check if using args works and then remove this line
    from Analyses.MetabolonAnalysis import *
    M = MetabolonEMPairwiseCorrelationWriter()
    if command_args.EMGenesPath is not None:
        M.A = Utils.Load(command_args.ppmetPath)
        M.B = Utils.Load(command_args.EMGenesPath)
    else:
        if command_args.data == 'KO':
            with open(command_args.ko_path, 'rb') as handle:
                M.B = pickle.load(handle)
            with open(command_args.ppmetPath, 'rb') as handle:
                M.A = pickle.load(handle)
        elif command_args.data == 'EMGenes':
            M.run()
        
    if command_args.output_dir is None and command_args.analysisType == 'single':
        command_args.output_dir = METABOLON_DIR + '/' + command_args.data + '_single_analysis/'
    elif command_args.output_dir is None and command_args.analysisType == 'shuffle':
        command_args.output_dir = METABOLON_DIR + '/' + command_args.data + '_shuffled_analysis/'
    
    if not os.path.exists(command_args.output_dir):
        os.makedirs(command_args.output_dir)
    
    with open(command_args.output_dir + '/args' + str(datetime.now()), 'w') as handle:
        for arg in vars(command_args):
            handle.write(str(arg) + '\t' + str(getattr(command_args, arg)) + '\n')
            
    if command_args.remove_nonKEGG_metabs:
        M.A = remove_non_kegg_metabolites(M.A, command_args.metabsPath)
        
    if command_args.remove_high_corr_metabs:
        M.A = remove_highly_correlated_metabolites(M.A)
    
    if command_args.from_file is not None:
        M.A = filter_metabolites_from_file(M.A, command_args.from_file)
        
    if command_args.analysisType == 'single':
        print (now() + " - uploading jobs to q...")
        with qp(jobname = 'BldSng', q=['himem7.q'], mem_def = command_args.memory_per_job, 
                trds_def = 1, max_u = 120, tryrerun = True, deleteCSHwithnoerr = True) as q:
            os.chdir(METABOLON_DIR)
            q.startpermanentrun()
            res = upload_jobs(q, M.A, M.B, 
                              command_args.output_dir, 
                              command_args.minimal_corr_pair, 
                              command_args.directional_pval, 
                              command_args.n_cols_per_job,
                              command_args.write_corr)
        concat_temp_corr_dfs(command_args.output_dir, command_args.n_cols_per_job, 
                             command_args.write_corr)
    elif command_args.analysisType == 'shuffle':
        with qp(jobname = 'BldShf', q=['himem7.q'], mem_def = command_args.memory_per_job, 
                trds_def = 1, max_u = 320, 
                tryrerun = True, deleteCSHwithnoerr = True) as q:
#             M.B = M.B.iloc[:,0:10000]
#             M.A = M.A.iloc[:,0:2]
            os.chdir(METABOLON_DIR)
            q.startpermanentrun()
            res = upload_jobs_shuffled(q, M.A, M.B, 
                                       command_args.output_dir, 
                                       command_args.minimal_corr_pair, 
                                       command_args.directional_pval,
                                       command_args.num_of_shuffles,
                                       command_args.write_corr,
                                       command_args.pass_path_to_q)

    return

    
    
        
if __name__ == "__main__":
    sethandlers()
    main()