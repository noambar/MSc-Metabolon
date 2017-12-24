###################################################################################################
# File: MetabolonGWAS.py
# Version: 0.0
# Date: 4.9.2017
# Noam Bar, noam.bar@weizmann.ac.il
#
# 
# Python version: 2.7
###################################################################################################
from Analyses.KEGGDataBuilder import *
import pandas as pd
from scipy.stats.stats import mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pickle
import os
from pandas import concat, read_csv, Series
from addloglevels import sethandlers
import argparse
from stats.HigherLevelRanksum import directed_mannwhitneyu
from queue.qp import qp
import re

# KO = False
# SHUFFELED_KO = False

METABOLON_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Metabolon/'
# WORKING_DIR = 'EMGenes_temp_corr_files/'
# DS_COLS = 'bac_genes'
KEGG_DB = ['reaction', 'pathway', 'module']

# if KO:
#     WORKING_DIR = '/ko_temp_corr_files/'
#     DS_COLS = 'kos'
# if SHUFFELED_KO and KO:
#     WORKING_DIR = '/ko_shuffeled_temp_corr_files/'
# WORKING_DIR = METABOLON_DIR + WORKING_DIR
# 
# EM_GENES_ALL_CORRS = WORKING_DIR + 'corrs.csv'
# EM_GENES_ALL_PVALS = WORKING_DIR + 'pvals.csv'

# METABOLON_FIG_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Metabolon/figs/'
# METABS_PATH = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Metabolon/tmp_files/metabs.df'

# GSEA_PATH = WORKING_DIR + 'GSEA.df'
# GSEA_SHUFFELED_PATH = WORKING_DIR + 'GSEA_shuffeled.df'

POSSIBLE_KEGG_DB = ['reaction', 'pathway', 'module', 'ko']
POSSIBLE_DOWNSTREAM_COL = ['bac_genes', 'ko']
POSSIBLE_STAT_TESTS = {'directed_mannwhitneyu':directed_mannwhitneyu, 'mannwhitneyu':mannwhitneyu}

# FILES_PER_BATCH = 100

class MetabolonKEGGGWAS():
    
    # pass as args:
    # working dir, ko, shuffled, ds_cols, pvals path, fig paht, metabs path, gsea path
     # add main with args
    def __init__(self, work_dir, ds_cols, gsea, metabsPath, stat_test, kos_per_compound, 
                 compounds_per_ko, plot = False, take_kdb_from_file = None,
                 take_metab_kdb_pair_from_file = None):
    #         self._EM_genes_all_corrs_path = EM_GENES_ALL_CORRS
    #         self._EM_genes_all_pvals_path = EM_GENES_ALL_PVALS
        self.kegg = KEGG_data_holder()
        self.kegg._merge_dicts('ko', 'compound')
        self.kegg._merge_dicts('compound', 'ko')
        self.work_dir = work_dir
        self.ds_cols = ds_cols
    #         self.kegg = kegg
        with open(metabsPath, 'rb') as handle:
            self.metabs = pickle.load(handle)
        self._shuffeled = False
        self._gsea_path = work_dir + '/' + gsea
        self._stat_test = POSSIBLE_STAT_TESTS[stat_test]
        self._stat_test_name = stat_test
        self._kos_per_compound = kos_per_compound
        self._compounds_per_ko = compounds_per_ko
        
        if plot:
            self.plot_some_data()
    #             self._get_corrs_matrix()
        if take_kdb_from_file is not None:
            with open(take_kdb_from_file, 'r') as handle:
                lines = handle.readlines()
            self.kdbs_from_file = [k.split('\n')[0] for k in lines]
        else:
            self.kdbs_from_file = None
            
        if take_metab_kdb_pair_from_file is not None:
            with open(take_metab_kdb_pair_from_file, 'r') as handle:
                lines = handle.readlines()
            lines = [line.split('\n')[0] for line in lines]
            self.metab_kdb_pair_from_file = {}
            for line in lines:
                fs = line.split('\t')
                if fs[0] in self.metab_kdb_pair_from_file:
                    self.metab_kdb_pair_from_file[fs[0]] += [fs[1]]
                else:
                    self.metab_kdb_pair_from_file[fs[0]] = [fs[1]]
        else:
            self.metab_kdb_pair_from_file = None
        
    
    def _get_pvals_matrix(self, pval_path):
        if not os.path.exists(pval_path):
            print 'no pval file'
            return
        print (str(datetime.now()) + " - [MetabolonKEGGGWAS._get_pvals_matrix] Read pvals from csv")
        self._EMGenes_pvals = read_csv(pval_path)
        print (str(datetime.now()) + " - [MetabolonKEGGGWAS._get_pvals_matrix] Done reading pvals from csv")
        self._EMGenes_pvals = self._EMGenes_pvals.set_index('Unnamed: 0').T #TODO: check if this helps
        self._EMGenes_pvals = concat((self._EMGenes_pvals.reset_index(drop = True), \
                                self._EMGenes_pvals.reset_index()['index'].apply(Series).\
                                    rename(columns = {0:'CompID'})), axis = 1)
        try:
            self._EMGenes_pvals['CompID'] = self._EMGenes_pvals.CompID.astype(int)
        except:
            # case the index is not of type int, remove it
            self._EMGenes_pvals = self._EMGenes_pvals.iloc[[int(float(x)) == float(x) for x in 
                                                            self._EMGenes_pvals.CompID.values],:]
            self._EMGenes_pvals['CompID'] = self._EMGenes_pvals.CompID.astype(int)
        self._EMGenes_pvals = self._EMGenes_pvals.set_index('CompID')
        print (str(datetime.now()) + " - [MetabolonKEGGGWAS._get_pvals_matrix] Done.")
        
    def _get_corrs_matrix(self, corr_path):
        print (str(datetime.now()) + " - [MetabolonKEGGGWAS._get_corrs_matrix] Read corrs from csv")
        self._EMGenes_corrs = read_csv(corr_path)
        print (str(datetime.now()) + " - [MetabolonKEGGGWAS._get_corrs_matrix] Done reading corrs from csv")
        self._EMGenes_corrs = self._EMGenes_corrs.set_index('Unnamed: 0').T #TODO: check if this helps
        self._EMGenes_corrs = concat((self._EMGenes_corrs.reset_index(drop = True), \
                                self._EMGenes_corrs.reset_index()['index'].apply(Series).\
                                    rename(columns = {0:'CompID'})), axis = 1)
        self._EMGenes_corrs['CompID'] = self._EMGenes_corrs.CompID.astype(int)
        self._EMGenes_corrs = self._EMGenes_corrs.set_index('CompID')
        print (str(datetime.now()) + " - [MetabolonKEGGGWAS._get_corrs_matrix] Done.")
        
        
    def metabolite_to_compound(self, metabolite):
        try:
            compound = str(self.metabs.loc[metabolite].KEGG)
            if compound == 'nan':
                return None
            return compound
        except:
            print ("metabolite " + metabolite + " does not exist in our dataset.")
            return None
    

    def plot_some_data(self):
        METABOLON_FIG_PATH = self.work_dir + '/fig/'
        def autolabel(rects, axx):
            """
            Attach a text label above each bar displaying its height
            """
            for rect in rects:
                height = rect.get_height()
                axx.text(rect.get_x() + rect.get_width()/2., 1.05*height, '%d' % int(height), 
                        ha='center', va='bottom')
            return axx
        
        def clip_X_percent(arr, X, top='top'):
            arr = arr[arr>0]
            if top == 'top':
                arr = arr.clip(max=np.percentile(arr, X))
            elif top == 'bottom':
                arr = arr.clip(min=np.percentile(arr, X))
            return arr
            
        # how many compounds with KEGG?
        sns.set_style("darkgrid")
        
        num_of_metabs = self.metabs.shape[0]
        only_kegg = self.metabs.KEGG[self.metabs.KEGG.isnull() == False]
        
        # out of the compounds with kegg, how many have some data? (reaction/pathway/module)
        # and many have downstream bacterial genes
        n_comp_w_db = {}
        n_comp_w_db_w_bg = {}
        db_colors = {'reaction': 'blue', 'pathway':'red', 'module':'green'}
        for db in ['reaction', 'pathway', 'module']:
            comp_w_db = [self.get_kegg_dbs_per_compound(i, db) for i in only_kegg.index]
            n_comp_w_db[db] = len([1 for i in comp_w_db if len(i) > 0])
            comp_w_db_w_bac_genes = np.array([sum([len(self.kegg.get_link_dicts(db + '_ko', 'ko_bacgene', r)) for r in c]) 
                                    for c in comp_w_db])

            n_comp_w_db_w_bg[db] = len([1 for i in comp_w_db_w_bac_genes if i > 0])

            # histograms of number of reactions/pathways/modules per compound
            n_db_per_compound = np.array([len(r) for r in comp_w_db if len(r) > 0])
            n_db_per_compound_95 = clip_X_percent(n_db_per_compound, 95, 'top')
            fig, ax = plt.subplots()
            ax.hist(n_db_per_compound_95, bins=15, color = db_colors[db])
            ax.set_ylabel('Number of compounds')
            ax.set_xlabel('Number of ' + db + 's')
            ax.set_title('Number of ' + db + 's per compound clipped at 95%')
            labels = ax.get_xticks().tolist()
            labels[-1] = '>' + str(labels[-1])
            ax.set_xticklabels(labels)
            plt.savefig(METABOLON_FIG_PATH + db + 's_per_compounds.png', bbox_inches='tight')
            
            # foreach kegg database, distribution of number of bac genes
            db_set = set([r for sublist in comp_w_db for r in sublist])
            bac_gene_per_db = np.array([len(self.kegg.get_link_dicts(db + '_ko', 'ko_bacgene', r)) for r in db_set])
            bac_gene_per_db_95 = clip_X_percent(bac_gene_per_db, 95, 'top')
            fig, ax = plt.subplots()
            ax.hist(bac_gene_per_db_95, bins=50, color = db_colors[db])
            ax.set_ylabel('Number of ' + db + 's')
            ax.set_xlabel('Number of bacterial genes')
            ax.set_title('Number of bacterial genes per ' + db + ' clipped at 95%')
            labels = ax.get_xticks().tolist()
            labels[-1] = '>' + str(labels[-1])
            ax.set_xticklabels(labels)
            plt.savefig(METABOLON_FIG_PATH + 'bac_gene_per_' + db + '.png', bbox_inches='tight')
            
            # KOs per kegg db
            ko_per_db = np.array([len(self.kegg.get_dicts()[db + '_ko'][d]) for d in self.kegg.get_dicts()[db + '_ko']])
            ko_per_db_95 = clip_X_percent(ko_per_db, 95, 'top')
            fig, ax = plt.subplots()
            ax.hist(ko_per_db_95, bins=15, color = db_colors[db])
            ax.set_ylabel('Number of ' + db + 's')
            ax.set_xlabel('Number of KOs')
            ax.set_title('Number of KOs per ' + db + ' clipped at 95%')
            labels = ax.get_xticks().tolist()
            labels[-1] = '>' + str(labels[-1])
            ax.set_xticklabels(labels)
            plt.savefig(METABOLON_FIG_PATH + 'KOs_per_' + db + '.png', bbox_inches='tight')
            
        # kegg data of compounds
        num_of_kegg_metabs = only_kegg.shape[0]
        num_of_kegg_compounds = len([c for c in only_kegg if c[0] == 'C'])
        num_of_kegg_drugs = len([d for d in only_kegg if d[0] == 'D'])
        nums = {'total metabs':num_of_metabs, 'metabs with KEGG id':num_of_kegg_metabs,
                'KEGG compounds':num_of_kegg_compounds, 'KEGG drugs':num_of_kegg_drugs}
        ind = np.arange(4)
        width = 0.7       
        fig, ax = plt.subplots()
        rects1 = ax.bar(ind, nums.values(), width, color='r')
        ax.set_ylabel('Number of metabolites')
        ax.set_title('Metabolites by KEGG')
        ax.set_xticks(ind + width / 2)
        ax.set_xticklabels(nums.keys())               
        ax = autolabel(rects1, ax)
        plt.savefig(METABOLON_FIG_PATH + 'metabolites_kegg_description.png', bbox_inches='tight')
        
        ind = np.arange(3)
        width = 0.35
        fig, ax = plt.subplots()
        rects1 = ax.bar(ind, n_comp_w_db.values(), width, color='r')
        rects2 = ax.bar(ind + width, n_comp_w_db_w_bg.values(), width, color='b')
        ax.set_ylabel('Number of metabolites')
        ax.set_title('Metabolites KEGG potential')
        ax.set_xticks(ind + width)
        ax.set_xticklabels(n_comp_w_db.keys())
        ax.legend((rects1[0], rects2[0]), ('compounds with KEGG data', 'compounds with bac genes'))              
        ax = autolabel(rects1, ax)
        ax = autolabel(rects2, ax)
        plt.savefig(METABOLON_FIG_PATH + 'metabolites_kegg_potential.png', bbox_inches='tight')
        
        # genes per KO
        bac_genes_per_ko = np.array([len(self.kegg.get_dicts()['ko_bacgene'][d]) for d in self.kegg.get_dicts()['ko_bacgene']])
        bac_genes_per_ko_95 = clip_X_percent(bac_genes_per_ko, 95, 'top')
        fig, ax = plt.subplots()
        ax.hist(bac_genes_per_ko_95, bins=50, color = 'grey')
        ax.set_ylabel('Number of KOs')
        ax.set_xlabel('Number of bacterial genes')
        ax.set_title('Number of bacterial genes per KO clipped at 95%')
        labels = ax.get_xticks().tolist()
        labels[-1] = '>' + str(labels[-1])
        ax.set_xticklabels(labels)
        plt.savefig(METABOLON_FIG_PATH + 'bacterial_genes_per_KO.png', bbox_inches='tight')
        
    def single_rank_sum_test(self, compound, genes, kegg, KEGG_database, gsea_id, min_dsc):
        # compute the rank sum test
        comp_pvals = self._EMGenes_pvals.loc[compound]
        try:
            db_pvals = comp_pvals[genes]
        except:
#             print (str(datetime.now()) + " - [MetabolonKEGGGWAS.single_rank_sum_test] No genes in index")
            return
        other_pvals = comp_pvals[self._EMGenes_pvals.columns.difference(genes)]
        db_pvals = db_pvals.dropna()
        other_pvals = other_pvals.dropna()
        if db_pvals.shape[0] < min_dsc or other_pvals.shape[0] == 0:
#             print (str(datetime.now()) + " - [MetabolonKEGGGWAS.single_rank_sum_test] There are not\
#              enough values in order to compute the test.")
            return
        res = self._stat_test(db_pvals.values, other_pvals.values)
        self.result_df.loc[self.result_df.shape[0]] = [compound, self.metabolite_to_compound(compound), 
                                                       KEGG_database, kegg, len(genes), db_pvals.shape[0], 
                                                       other_pvals.shape[0], res[1], res[0], 
                                                       self._stat_test_name, self.ds_cols, gsea_id]
        return
        
    def get_kegg_dbs_per_compound(self, compound, KEGG_database):
        
        kegg_compound = self.metabolite_to_compound(compound)
        if (kegg_compound is None):
#             print (str(datetime.now()) + " - following compound has no KEGG identifier: "  + str(compound))
            return []
        #drugs KEGG identifier start with 'D', compounds start with 'C'
        if kegg_compound[0] == 'D':
            c_d = 'drug'
        else:
            c_d = 'compound'
        dic_name = c_d + '_' + KEGG_database
        kegg_pref = 'cpd:' if c_d == 'compound' else 'dr:'
        kegg_compound = kegg_pref + kegg_compound
        if dic_name not in self.kegg.get_dicts():
            print (str(datetime.now()) + " - The following link is unavailable: "  + dic_name)
            return []
        if kegg_compound not in self.kegg.get_dicts()[dic_name]:
#             print (str(datetime.now()) + " - " + kegg_compound + " not in "  + dic_name)
            return []
        return self.kegg.get_dicts()[dic_name][kegg_compound]
        
    def rank_sum_test_per_compound(self, compound, KEGG_database, gsea_id, 
                                   min_dsc, max_dsc, max_kos_per_db):
        # for each kegg database call single_rank_sum_test
        if self.metab_kdb_pair_from_file is not None:
            if str(compound) not in self.metab_kdb_pair_from_file:
                return
            keggs = self.metab_kdb_pair_from_file[str(compound)]
        elif (KEGG_database == 'ko'):
            # screening metabolites with too many kos
            kegg_compound = 'cpd:' + self.metabolite_to_compound(compound)
            if kegg_compound not in self.kegg.get_dicts()['compound_ko']:
                return
            if len(self.kegg.get_dicts()['compound_ko'][kegg_compound]) > self._kos_per_compound:
                print compound
                return
            keggs = []
            for kdb in KEGG_DB:
                kdbs = self.get_kegg_dbs_per_compound(compound, kdb)
                for k in kdbs:
                    if k in self.kegg.get_dicts()[kdb + '_ko']:
                        kos = self.kegg.get_dicts()[kdb + '_ko'][k]
                        if len(kos) <= max_kos_per_db:
                            keggs += kos
            keggs = list(set(keggs))

        else:        
            keggs = self.get_kegg_dbs_per_compound(compound, KEGG_database)
            
#         print (str(datetime.now()) + " - [MetabolonKEGGGWAS.rank_sum_test_per_compound] kegg dbs: " + str(len(keggs)))
        for k in keggs:
            if self.kdbs_from_file is not None:
                if k not in self.kdbs_from_file:
                    continue

            if self.ds_cols == 'bac_genes':
                if (KEGG_database == 'ko'):
                    # check that KO doesn't have too many compounds
                    if len(self.kegg.get_dicts()['ko_compound'][k]) > self._compounds_per_ko:
                        print k
                        continue
                    
                    if k in self.kegg.get_dicts()['ko_bacgene']:
                        genes = self.kegg.get_dicts()['ko_bacgene'][k]
                    else:
                        genes = []
                else:
                    genes = self.kegg.get_link_dicts(KEGG_database + '_ko', 'ko_bacgene', k)
            elif self.ds_cols == 'kos':
                if k in self.kegg.get_dicts()[KEGG_database + '_ko']:
                    genes = self.kegg.get_dicts()[KEGG_database + '_ko'][k]
                    genes = [ko.split('ko:')[1] for ko in genes]
                else:
                    genes = []
            if len(genes) >= min_dsc and len(genes) <= max_dsc:
                self.single_rank_sum_test(compound, genes, k, KEGG_database, gsea_id, min_dsc)
            else:
                continue
        return
    
    def all_rank_sum_tests(self, KEGG_database, gsea_id, min_dsc, max_dsc, max_kos_per_db):
        # call rank_sum_test_per_compound for each compound
        print (str(datetime.now()) + " - [MetabolonKEGGGWAS.all_rank_sum_tests] Start: " + str(gsea_id))
        if os.path.exists(self._gsea_path):
            with open(self._gsea_path, 'rb') as handle:
                self.result_df = pickle.load(handle)
        else:
            self.result_df = pd.DataFrame(columns=['CompID', 'KEGGID', 'KEGG_db_type', 'KEGG_db_ID', 
                                               'n_genes_in_kegg', 'n_genes_in_test', 'n_genes_rest', 
                                               'test_pval', 'test_statistic', 'stat_test', 'downstream_cols', 'GSEA_ID'])
        comp_counter = self._EMGenes_pvals.shape[0]
        for comp in self._EMGenes_pvals.index:
            comp_counter += (-1)
            if comp_counter % 10 == 0:   
                print (str(datetime.now()) + " - [MetabolonKEGGGWAS.all_rank_sum_tests] compID: " + 
                      str(comp) + ", " + str(comp_counter) + " to go...")
            self.rank_sum_test_per_compound(int(comp), KEGG_database, gsea_id, min_dsc, max_dsc, max_kos_per_db)
            
        with open(self._gsea_path, 'wb') as handle:
            pickle.dump(self.result_df, handle, protocol=pickle.HIGHEST_PROTOCOL)
        print (str(datetime.now()) + " - [MetabolonKEGGGWAS.all_rank_sum_tests] End")
        
        
    def perform_shuffeled_rank_sum_test(self, KEGG_database, downstream_cols, back_to_normal = False):
        print (str(datetime.now()) + " - [MetabolonKEGGGWAS.perform_shuffeled_rank_sum_test] Start")
        if (self._shuffeled == False):
            self._orig_gsea_path = self._gsea_path
            self._gsea_path = self._gsea_path + 'suffle'
            self._orig_cols = self._EMGenes_pvals.columns
            print (str(datetime.now()) + " - [MetabolonKEGGGWAS.perform_shuffeled_rank_sum_test] Permuting columns")
            new_cols = np.random.permutation(self._EMGenes_pvals.columns)
            self._EMGenes_pvals.columns = new_cols
            self._shuffeled = True
        self.all_rank_sum_tests(KEGG_database)
        if (back_to_normal):
            self._EMGenes_pvals.columns = self._orig_cols
            self._gsea_path = self._orig_gsea_path
            self._shuffeled = False
        return
    
    
    
def divide_jobs(command_args, files_list, idx):
    mkg = MetabolonKEGGGWAS(command_args.work_dir, command_args.down_stream_cols, command_args.gsea, 
                            command_args.metabsPath, command_args.stat_test,
                            command_args.kos_per_compound, command_args.compounds_per_ko,
                            command_args.plot, command_args.take_kdb_from_file,
                            command_args.take_metab_kdb_pair_from_file)
    mkg._gsea_path = mkg._gsea_path + '_' + str(idx)
    for f in files_list:
        mkg._get_pvals_matrix(command_args.work_dir + '/' + f)
        if command_args.keggdb == 'all':
            for kdb in POSSIBLE_KEGG_DB:
                mkg.all_rank_sum_tests(kdb, f, command_args.minimal_downstream_cols, 
                                       command_args.maximal_downstream_cols, 
                                       command_args.maximal_kos_per_db)
        else:
            mkg.all_rank_sum_tests(command_args.keggdb, f, command_args.minimal_downstream_cols, 
                                    command_args.maximal_downstream_cols, 
                                    command_args.maximal_kos_per_db)
    return



def concat_gsea_files(work_dir, prefix):
    print(str(datetime.now()) + " - [concat_gsea_files] start")
    gsea_files = os.listdir(work_dir)
    gsea_files = [f for f in gsea_files if len(f) >= len(prefix) and f[:len(prefix)] == prefix]
    final_df = pd.DataFrame()
    for f in gsea_files:
        print f
        f_name = work_dir + '/' + f
        temp_df = Utils.Load(f_name)
        final_df = concat((final_df, temp_df), axis=0)
        os.remove(f_name)
    Utils.Write(work_dir + '/' + prefix, final_df)
    print(str(datetime.now()) + " - [concat_gsea_files] end")
    return

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    '''
    return [atoi(c) for c in re.split('(\d+)', text)]
        


def main():
    """
    Compute all correlations between EM genes matrix and metabolites and write to disk.
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument('work_dir', help='Path to working directory', type=str, default=None)
    parser.add_argument('-analysisType', help='What type of analysis to perform: normal or shuffle', type=str, default='normal')
    parser.add_argument('-keggdb', help='Which KEGG DB to consider: reaction, pathway, module, ko, all', type=str, default='ko')
    parser.add_argument('-plot', help='Whether to call the plotting function', type=bool, default=False)
    parser.add_argument('-gsea', help='Name of GSEA results file in working directory', type=str, default='GSEA')
    parser.add_argument('-metabsPath', help='Path to metabolites file', type=str, default=METABOLON_DIR + '/tmp_files/metabs.df')
    parser.add_argument('-stat_test', help='Type of statistical test: mannwhitneyu, directed_mannwhitneyu', type=str, default='directed_mannwhitneyu')
    parser.add_argument('-dsc', '--down_stream_cols', help='What type of data: bac_genes or ko', type=str, default='bac_genes')
    parser.add_argument('-tkdbff', '--take_kdb_from_file', help='Path to KEGG DB', type=str, default=None)
    parser.add_argument('-tkprff', '--take_metab_kdb_pair_from_file', help='Path to metabolite-KEGGDB pairs file', type=str, default=None)
    parser.add_argument('-min_dsc', '--minimal_downstream_cols', help='Minimal number downstream columns in test', type=int, default=5)
    parser.add_argument('-max_dsc', '--maximal_downstream_cols', help='Maximal number downstream columns in test', type=int, default=10000)
    parser.add_argument('-max_kos_db', '--maximal_kos_per_db', help='Maximal number KOs per KEGG DB', type=int, default=100)
    parser.add_argument('-fpb', '--files_per_batch', help='Number of files to parse in each job', type=int, default=100)
    parser.add_argument('-kpc', '--kos_per_compound', help='Maximal number of KOs per compound', type=int, default=2000)
    parser.add_argument('-cpk', '--compounds_per_ko', help='Maximal number of compound per KO', type=int, default=500)
    
    command_args = parser.parse_args()

    if command_args.work_dir is None:
        print(str(datetime.now()) + " - [main] work_dir must be passed.")
        return
    if command_args.analysisType != 'normal' and command_args.analysisType != 'shuffle':
        print(str(datetime.now()) + " - [main] analysisType must be normal or shuffle.")
        return
    if command_args.down_stream_cols not in POSSIBLE_DOWNSTREAM_COL:
        print(str(datetime.now()) + " - [main] down_stream_cols must be bac_genes or ko.")
        return
    if command_args.keggdb != 'all' and command_args.keggdb not in POSSIBLE_KEGG_DB:
        print(str(datetime.now()) + " - [main] Illegal KEGG DB.")
        return
    if command_args.stat_test not in POSSIBLE_STAT_TESTS:
        print(str(datetime.now()) + " - [main] Illegal stat test.")
        return
    
    with open(command_args.work_dir + '/args' + str(datetime.now()), 'w') as handle:
        for arg in vars(command_args):
            handle.write(str(arg) + '\t' + str(getattr(command_args, arg)) + '\n')
            
#     mkg = MetabolonKEGGGWAS(command_args.work_dir, command_args.down_stream_cols, command_args.gsea, 
#                             command_args.metabsPath, command_args.stat_test, 
#                             command_args.plot, command_args.take_kdb_from_file,
#                             command_args.take_metab_kdb_pair_from_file)
    
    if command_args.analysisType == 'normal':
        pval_files = os.listdir(command_args.work_dir)
        pval_files = [f for f in pval_files if len(f) > 4 and f[0:5] == 'pvals']
        pval_files.sort(key=natural_keys)
        if len(pval_files) > command_args.files_per_batch:
            with qp('MetGWAS', delay_sec=1, mem_def='10G', q = ['himem7.q'], trds_def = 1, 
                    max_u = 120, tryrerun=True) as q:
                os.chdir(METABOLON_DIR)
                waiton = []
                q.startpermanentrun()
                for i in range(0, len(pval_files), command_args.files_per_batch):
                    print(str(datetime.now()) + " - [main] uploading job " + str(i) + " to q")
                    waiton.append(q.method(divide_jobs, 
                       (command_args, pval_files[i:i+command_args.files_per_batch], i)))
                res = q.waitforresults(waiton)
            print(str(datetime.now()) + " - [main] Results are back")
            concat_gsea_files(command_args.work_dir, command_args.gsea)
            return

        else:
            mkg = MetabolonKEGGGWAS(command_args.work_dir, command_args.down_stream_cols, command_args.gsea, 
                                    command_args.metabsPath, command_args.stat_test,
                                    command_args.kos_per_compound, command_args.compounds_per_ko,
                                    command_args.plot, command_args.take_kdb_from_file,
                                    command_args.take_metab_kdb_pair_from_file)
            for f in pval_files:
                mkg._get_pvals_matrix(command_args.work_dir + '/' + f)
                if command_args.keggdb == 'all':
                    for kdb in POSSIBLE_KEGG_DB:
                        mkg.all_rank_sum_tests(kdb, f, command_args.minimal_downstream_cols, 
                                               command_args.maximal_downstream_cols, 
                                               command_args.maximal_kos_per_db)
                else:
                    mkg.all_rank_sum_tests(command_args.keggdb, f, command_args.minimal_downstream_cols, 
                                           command_args.maximal_downstream_cols, 
                                           command_args.maximal_kos_per_db)
    elif command_args.analysisType == 'shuffle':
        # no much reason to use this, since shuffling is done before hand
        pass




# mkg = MetabolonKEGGGWAS()
# for kdb in KEGG_DB:
#     mkg.perform_shuffeled_rank_sum_test(kdb, DS_COLS, 0)
# mkg.all_rank_sum_tests('ko')
# mkg.all_rank_sum_tests('reaction', downstream_cols=DS_COLS)
# mkg.all_rank_sum_tests('pathway', downstream_cols=DS_COLS)
# mkg.all_rank_sum_tests('module', downstream_cols=DS_COLS)
# mkg.perform_shuffeled_rank_sum_test('reaction')
# mkg.perform_shuffeled_rank_sum_test('pathway')
# mkg.perform_shuffeled_rank_sum_test('module')

# running shuffeled KOs (turn SHUFFELED_KO and KO to True)
# mkg = MetabolonKEGGGWAS()
# pval_files = os.listdir(WORKING_DIR)
# pval_files = [f for f in pval_files if f[0] == 'p']
# for f in pval_files:
#     mkg._get_pvals_matrix(WORKING_DIR + f)
#     for kdb in KEGG_DB:
#         mkg.all_rank_sum_tests(kdb, DS_COLS, f)
    
    
if __name__ == "__main__":
    sethandlers()
    main()