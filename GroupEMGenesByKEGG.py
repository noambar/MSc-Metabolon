###################################################################################################
# File: GroupEMGenesByKO.py
# Version: 0.1
# Date: 3.10.2017
# Noam Bar, noam.bar@weizmann.ac.il
#
# 
# Python version: 2.7
###################################################################################################

from Analyses.KEGGDataBuilder import *
# from scipy.stats.stats import mannwhitneyu, spearmanr, pearsonr
# from datetime import datetime
# import pickle
# import numpy as np
import pandas as pd
import numpy as np
from pandas import concat
from Analyses.grouping import semi_exahustive_strict_cor_grouping, newalg
from addloglevels import sethandlers
from queue.qp import qp
import os
import argparse

METABOLON_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Metabolon/'
POSSIBLE_KEGG_DB = ['reaction', 'pathway', 'module', 'ko']
POSSIBLE_METHODS = ['newalg', 'sescg']
REASONABLE_CORRELATIONS = ['spearman', 'pearson']
REASONABLE_GROUPING_TYPES = ['sum', 'rep']
                
def _grouping_with_sescg(mat, k, work_dir, corr_thresh, abscorr, presence_thresh, 
                         corr_presence_thresh, first_split_n, split_n, corrmethod, grouping_type):
    print(str(datetime.now()) + " - [_grouping_with_sescg] start - " + k)
    with qp('grping', delay_sec=1, tryrerun=True, mem_def='2G', trds_def=2, q = ['himem7.q']) as q:
        q.startpermanentrun()
        ex_grps = semi_exahustive_strict_cor_grouping(mat, corr_thresh, abscorr, presence_thresh,
                                                       corr_presence_thresh, first_split_n,
                                                       split_n, corrmethod, q)
    cormet = mat.corr(corrmethod, corr_presence_thresh)
    def select_rep(grp):
        if len(grp) == 1:
            return grp[0]
        grp_nnl = mat[grp].notnull().sum()
        return cormet.loc[grp_nnl[grp_nnl == grp_nnl.max()].index, 
                          grp_nnl[grp_nnl == grp_nnl.max()].index].sum(1).argmax()
    met_lg = concat((mat[select_rep(grp)] for grp in ex_grps), axis = 1)
    Utils.Write(work_dir + '/' + k, met_lg)
    print(str(datetime.now()) + " - [_grouping_with_sescg] end - " + k)
    return met_lg

def _grouping_with_newalg(mat, k, work_dir, corr_thresh, abscorr, corr_presence_thresh, corrmethod, 
                          mirror, grouping_type):
    print(str(datetime.now()) + " - [_grouping_with_newalg] start - " + k)
    cormet = mat.corr(corrmethod, corr_presence_thresh)
    if abscorr:
        cormet = cormet.abs()
    # case where all genes are in correlation, the newalg fails
    if (cormet.values.ravel() > corr_thresh).sum() == cormet.shape[0] * cormet.shape[1]:
        ex_grps = [list(cormet.columns)]
    else:
        ex_grps = newalg(cormet, corr_thresh, mirror)
    def select_rep(grp):
        if len(grp) == 1:
            return grp[0]
        grp_nnl = mat[grp].notnull().sum()
        return cormet.loc[grp_nnl[grp_nnl == grp_nnl.max()].index, 
                          grp_nnl[grp_nnl == grp_nnl.max()].index].sum(1).argmax()

    if grouping_type == 'rep':
        met_lg = concat((mat[select_rep(grp)] for grp in ex_grps), axis = 1)
    elif grouping_type == 'sum':
        mat = mat.apply(lambda x: 10**x)
        mat = mat.fillna(0)
        met_lg = concat((pd.DataFrame(mat[grp].mean(1).values, index=mat.index, columns=[select_rep(grp)]) for grp in ex_grps), axis = 1)
        met_lg[met_lg == 0] = np.nan
        met_lg = met_lg.apply(np.log10)
    Utils.Write(work_dir + '/' + k, met_lg)
    print(str(datetime.now()) + " - [_grouping_with_newalg] end - " + k)
    return met_lg

def _concat_output(work_dir, output_name, prefix, delete=False):
    print(str(datetime.now()) + " - [_concat_output] start")
    df_list = os.listdir(work_dir)
    df_list = [df for df in df_list if len(df) >= len(prefix) and df[:len(prefix)] == prefix]
    final_df = pd.DataFrame()
    i = 0
    b = len(df_list)/20.
    for df_path in df_list:
        with open(work_dir + '/' + df_path, 'rb') as handle:
            df = Utils.Load(handle)
        final_df = concat((final_df, df), axis=1)
        if i % b < 1:
            print str(int(100.*i/len(df_list))) + '%'
        i += 1
        if delete:
            os.remove(work_dir + '/' + df_path)
    Utils.Write(work_dir + '/' + output_name, final_df)
    print(str(datetime.now()) + " - [_concat_output] end")
    return
        
def main():
    print(str(datetime.now()) + " - [GroupEMGenesByKEGG:main] start")
    parser = argparse.ArgumentParser()
    parser.add_argument('-work_dir', help='Path to working directory', type=str, default=None)
    parser.add_argument('-output_name', help='Name of output dataframe file', type=str, default='EMGenes_group.df')
    parser.add_argument('-keggDB', help='Type of KEGG database to group by', type=str, default='ko')
    parser.add_argument('-min_size', help='Minimum size of KEGG DB to perform grouping on', type=int, default=20)
    parser.add_argument('-max_size', help='Minimum size of KEGG DB to perform grouping on', type=int, default=100000)
    parser.add_argument('-method', help='Which algorithm to use for grouping: newalg or sescg', type=str, default='newalg')
    parser.add_argument('-grouping_type', help='Whether to choose a representetive or to sum: sum or rep', type=str, default='rep')
    parser.add_argument('-corr_thresh', help='Correlation threshold to group by', type=float, default=0.75)
    parser.add_argument('-abscorr', help='Use absolute value of correlation', type=bool, default=True)
    parser.add_argument('-presence_thresh', help='Demand at least this amount of samples', type=int, default=100)
    parser.add_argument('-corr_presence_thresh', help='Threshold for computing correlation', type=int, default=50)
    parser.add_argument('-first_split_n', help='first split n', type=int, default=30)
    parser.add_argument('-split_n', help='split n', type=int, default=1)
    parser.add_argument('-corrmethod', help='Which correlation method to use', type=str, default='spearman')
    parser.add_argument('-mirror', help='mirror option in newalg', type=bool, default=False)
    parser.add_argument('-just_concat', help='Whether to just read and concat the output', type=bool, default=False)
    command_args = parser.parse_args()
    
    if command_args.work_dir is None:
        command_args.work_dir = METABOLON_DIR + '/Grouping_by_' + command_args.keggDB
    if not os.path.exists(command_args.work_dir):
        os.makedirs(command_args.work_dir)

    if command_args.keggDB not in POSSIBLE_KEGG_DB:
        print(str(datetime.now()) + " - [GroupEMGenesByKEGG:main] keggDB not legal")
        return
    
    if str(command_args.min_size).isdigit() is False or str(command_args.max_size).isdigit() is False:
        print(str(datetime.now()) + " - [GroupEMGenesByKEGG:main] min_size and max_size must be integers")
        return
    
    if command_args.method not in POSSIBLE_METHODS:
        print("method must be one of: " + ', '.join(POSSIBLE_METHODS))
        return
    
    if command_args.corr_thresh > 1 or command_args.corr_thresh < -1:
        print(str(datetime.now()) + " - [GroupEMGenesByKEGG:main] corr_thresh must be between -1 and 1")
        return
    
    if command_args.corrmethod not in REASONABLE_CORRELATIONS:
        print("corrmethod must be one of: " + ', '.join(REASONABLE_CORRELATIONS))
        return
        
    if command_args.grouping_type not in REASONABLE_GROUPING_TYPES:
        print("corrmethod must be one of: " + ', '.join(REASONABLE_GROUPING_TYPES))
        return
    
    if command_args.just_concat:
        _concat_output(command_args.work_dir, command_args.output_name, command_args.keggDB, False)
        return
    
    with open(command_args.work_dir + '/args' + str(datetime.now()), 'w') as handle:
        for arg in vars(command_args):
            handle.write(str(arg) + '\t' + str(getattr(command_args, arg)) + '\n')
    
    print(str(datetime.now()) + " - [GroupEMGenesByKEGG:main] build KEGG data holder")
    kegg = KEGG_data_holder()
    print(str(datetime.now()) + " - [GroupEMGenesByKEGG:main] Load EMGenes dataframe")
    from Analyses.MetabolonAnalysis import *
    M = MetabolonEMPairwiseCorrelationWriter()
    M.run()
    EMGenes = M.B
    
    print(str(datetime.now()) + " - [GroupEMGenesByKEGG:main] Divide EMGenes by " + command_args.keggDB)
    keggDict = {}
    EMGenes_columns_set = set(EMGenes.columns)
    if (command_args.keggDB == 'ko'):
        for k in kegg.get_dicts()['ko_bacgene'].keys():
            temp_genes = set(kegg.get_dicts()['ko_bacgene'][k])
            valid_genes = list(temp_genes.intersection(EMGenes_columns_set))
            if len(valid_genes) > 0:
                keggDict[k] = EMGenes[valid_genes]
                
                
    print(str(datetime.now()) + " - [GroupEMGenesByKEGG:main] uploading jobs to q...")
    
    with qp('m_grping', delay_sec=1, mem_def='1G', q = ['himem7.q'], trds_def = 1, max_u = 420, 
            tryrerun=True) as q_m:
        unused_genes = set()
        os.chdir(METABOLON_DIR)
        waiton = []
        q_m.startpermanentrun()
        for k in keggDict.keys():
            if keggDict[k].shape[1] < command_args.min_size or keggDict[k].shape[1] > command_args.max_size:
                unused_genes.update(keggDict[k].columns)
                continue
            if os.path.exists(command_args.work_dir + '/' + k):
                continue
            if keggDict[k].shape[1] > command_args.min_size:
                print k + ' - ' + str(keggDict[k].shape[1])
            if command_args.method == 'newalg':
                waiton.append(q_m.method(_grouping_with_newalg, 
                       (keggDict[k], k, command_args.work_dir, command_args.corr_thresh, 
                        command_args.abscorr, command_args.corr_presence_thresh, 
                        command_args.corrmethod, command_args.mirror, command_args.grouping_type)))
            elif command_args.method == 'sescg':
                waiton.append(q_m.method(_grouping_with_sescg, 
                       (keggDict[k], k, command_args.work_dir, command_args.corr_thresh, 
                        command_args.abscorr, command_args.presence_thresh, 
                        command_args.corr_presence_thresh, command_args.first_split_n, 
                        command_args.split_n, command_args.corrmethod, command_args.grouping_type)))
        res = q_m.waitforresults(waiton)
    unused_genes = list(unused_genes)
    Utils.Write(command_args.work_dir + '/ko:rest.df', EMGenes[unused_genes])
    
    print(str(datetime.now()) + " - [GroupEMGenesByKEGG:main] concatinating dataframes into one")
    _concat_output(command_args.work_dir, command_args.output_name, command_args.keggDB, True)
    print(str(datetime.now()) + " - [GroupEMGenesByKEGG:main] end")
    return
    
if __name__ == "__main__":
    sethandlers()
    main()