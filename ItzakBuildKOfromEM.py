###################################################################################################
# File: MetabolonAnalysis.py
# Version: 0.0
# Date: 14.8.2017
# Noam Bar, noam.bar@weizmann.ac.il
#
# 
# Python version: 2.7
###################################################################################################

from pandas import concat, read_csv, Series#, read_sql, read_pickle
from GPBackend.SingleExplorer import SingleExplorer
# from load_data import get_pnp_ids, get_preprocessed
from GPBackend.DataBuild.FeatureExposers.PisAndPTsExposer import FDCIDIndexedExposer
from abc import ABCMeta
from GPBackend.DataBuild.FeatureExposer import LoadAction#, register_features
from GPBackend.DataBuild.FeatureGroups.DynamicFeatureGroup import DynamicFeatureGroup
from GPBackend.ModDataBuild import SubDataBuilder, DataBuilder
from GPBackend.DataBuild.Transformers.FDSelector import FDSelect,\
    GenotekHandling, StudyTypeHandling, RemainingHandling
from GPBackend.DataBuild.FeatureExposers.MetagenomeMDExposer import MetagenomeMDExposer
from GPBackend.DataBuild.Prefilters.MinMetagenomeReadsFilter import MinReadFilter
# from GPBackend.DataBuild.Transformers.FillNA import FillNATransformer
# from GPBackend.DataBuild.Transformers.CapFromBelow import CapFromBelow
# from GPBackend.DataBuild.Transformers.TakeLog import TakeLog
# from GPBackend.DataBuild.Transformers.PCATransformer import PCATransformer
# from GPBackend.PipelineManipulate.TransformEstimators.GroupFeatureEstimators.TransformedFeatures import TransformedFeatures
from GPBackend.DataBuild.Transformer import Transformer
import numpy as np
from Caching.Caching import Signer
# from stats.statutils import intersected_v_method_mincount
# from scipy.stats.stats import spearmanr
# from os.path import join
# import matplotlib.pyplot as plt
# import seaborn as sns
import pandas as pd
import Utils

from datetime import datetime
import os
import pickle
from Analyses.grouping import semi_exahustive_strict_cor_grouping, strict_cor_grouping

METABOLON_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/'

EM_GENES_PATH_10 = METABOLON_DIR + 'EM_Genes_raw_first10.csv'
EM_GENES_PATH_RAW = '/net/mraid08/export/jafar/Microbiome/Analyses/Metabolon2/DFOut/EM_Genes_raw.csv'
# EM_GENES_PATH_RAW = EM_GENES_PATH_10
# AFTER_FIX_STUFF = METABOLON_DIR + 'EMGenes_after_Fix_stuff.csv'
# EM_GENES_PATH_READY = METABOLON_DIR + 'EMGenes_ready_483x3809385.csv'
EM_GENES_PATH_READY = METABOLON_DIR + 'EMGenes_ready_10PercentScreened.csv'
KEGG_COG_DICT_PATH = '/net/mraid08/export/jafar/Microbiome/Data/Databases/GeneSet2014/Annotation/KEGG_COG_dict.dat'
KO_EMGENES = METABOLON_DIR + 'KO_EMGenes_all.df'

# KEGG_COG_DICT_PATH = '/net/mraid08/export/jafar/Microbiome/Data/Databases/GeneSet2014/Annotation/KEGG_COG_dict.dat'


# Don't try to load the raw EM genes matrix if it is saved to memory
# if os.path.exists(EM_GENES_PATH_READY):
#     EM_GENES_PATH_RAW = EM_GENES_PATH_10

NA_FILLER = 1e-7

class AdHocBuilder(Signer):
    code_version = 2
    _composingitems = ()
    def __init__(self, df):
        self.df = df
    def getX(self):
        return self.df
    def drop_features(self, a, b):
        return self.df

class EMGenes_Exposer(FDCIDIndexedExposer):
    __metaclass__ = ABCMeta
    code_version = 2
    _composingitems = ()
    DataFramePath = EM_GENES_PATH_RAW
    loadmethod = LoadAction.LoadCsv
    G_All = DynamicFeatureGroup()
    def _postload(self):
        self._df = self._df.set_index('Unnamed: 0').T
        super(EMGenes_Exposer, self)._postload()
# register_features(EMGenes_Exposer)

class EMGenesConcat(Transformer):
    code_version = 1
    _composingitems = ()
    
    def transform(self, df):
        print (str(datetime.now()) + " - EMGenesConcat")
        print (str(datetime.now()) + " - Reading EM Genes matrix from csv")
        x = read_csv(EM_GENES_PATH_RAW)
        print (str(datetime.now()) + " - Done reading EM Genes matrix from csv")
        # assign genes as index and transpose
        x = x.set_index('Unnamed: 0').T
#         x = x.reindex(x.index.rename('index'))
        return x
#         print (str(datetime.now()) + " - Initial shape of EM genes matrix: " + str(x.shape))
#         # break index into FD and connectionID
#         x = concat((x.reset_index(drop = True), \
#                     x.reset_index()['index'].str.split('_').apply(Series).\
#                         rename(columns = {0:'FD', 1:'ConnectionID'})), axis = 1)
#         # remove samples were FD is null
#         x = x.loc[x['ConnectionID'] != 'nan']
#         
# #         x['ConnectionID'] = x.ConnectionID.astype(int)
#         x = x.sort_values(['ConnectionID', 'FD']).set_index(['ConnectionID', 'FD'])
#         # add some information regarding the samples
#         df2 = concat((df, x), axis = 1)
#         print (str(datetime.now()) + " - EM shape: " + str(df2.shape))
#         return df2
    
    
class FixStuff(Transformer):
    code_version = 1
    _composingitems = ()
    
    def transform(self, df):
#         if not os.path.exists(EM_GENES_PATH_READY):
        print (str(datetime.now()) + " - Fix stuff transform")
        print (str(datetime.now()) + " - Initial shape: " + str(df.shape))
        # remove columns added in EMGenesConcat
        df = df.iloc[:, 6:].dropna(1, how = 'all')
        print (str(datetime.now()) + " - After removing NA's: " + str(df.shape))
#         with open('/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Metabolon/tmp_files/pnp_ids.df', 'rb') as handle:
#             pnp_ids = pickle.load(handle)
#             
#         c = concat((pnp_ids[pnp_ids.StudyTypeID != 50],
#             pnp_ids[pnp_ids.StudyTypeID == 50].reset_index().groupby('ConnectionID').\
#                                                 apply(lambda x: x.loc[x.StorageDT.argmin()]).\
#                                                 set_index('CLIENT_IDENTIFIER')))
#         df = c[['ConnectionID']].merge(df, left_on = 'ConnectionID', right_index = True).drop('ConnectionID', axis = 1)
        
#         df = df.apply(np.log10)

        print (str(datetime.now()) + " - EM genes shape: " + str(df.shape))
        return df
    
class RemoveLowAbundantGenesTransformer(Transformer):
    code_version = 1
    _composingitems = ()
    
    def transform(self, df):    
        # taking only abundant genes
        print (str(datetime.now()) + " - RemoveLowAbundantGenesTransformer")    
        ten_per = df.shape[0] - (df.shape[0] / 10.)
        above10 = (df.isnull().sum() < ten_per).values
        df = df.iloc[:, above10]
        print (str(datetime.now()) + " - EM shape: " + str(df.shape))
        return df  
    
class BuildKOdfTransformer(Transformer):
    code_version = 1
    _composingitems = ()
    
    def transform(self, df):
        print (str(datetime.now()) + " - BuildKOdfTransformer")    
        if not os.path.exists(KO_EMGENES):
            print (str(datetime.now()) + " - EM genes: Fill na with 0 ...")
            df = df.fillna(0)
            KO_df = pd.DataFrame(index=df.index)
            print (str(datetime.now()) + " - Loading KEGG COG...")
            kegg = Utils.Load(KEGG_COG_DICT_PATH)
            print (str(datetime.now()) + " - Creating KEGG matrix...")
            for gene in df.columns:
                if gene in kegg:
                    kos = kegg[gene]['KEGG']
                    if kos != 'unknown':
                        kos = kos.split(';')
                        for ko in kos:
                            if ko in KO_df.columns:
                                KO_df[ko] += df[gene]
                            else:
                                KO_df[ko] = df[gene]
            print "KO matrix shape: " + str(KO_df.shape)
            KO_df[KO_df == 0] = np.nan
#             KO_df = KO_df.apply(np.log10)
            Utils.Write(KO_EMGENES, KO_df)                   
            df[df == 0] = np.nan
        print (str(datetime.now()) + " - EM shape: " + str(df.shape))
        assert 1==0
        return df 
    
class UniteCorrelatedEMGenesWithinKOTransformer(Transformer):
    code_version = 1
    _composingitems = ()
    
    def transform(self, df):
        print (str(datetime.now()) + " - UniteCorrelatedEMGenesWithinKOTransformer")        
        # foreach KO, go over genes in it, and compute correlations between them
        
        # Find a way to unite the genes smartly, maybe with hirrercal clustering
        print (str(datetime.now()) + " - writing ready EM genes to csv...")
        df.T.to_csv(EM_GENES_PATH_READY)
        print (str(datetime.now()) + " - EM shape: " + str(df.shape))
        return df
    
class RemoveGenesWithMoreThanOneKOTransfomer(Transformer):
    code_version = 1
    _composingitems = ()
    
    def transform(self, df):
        # remove bac genes which map to more than one KOs
        print (str(datetime.now()) + " - RemoveGenesWithMoreThanOneKOTransfomer")
        valid_genes = []
        kegg_cog = Utils.Load(KEGG_COG_DICT_PATH)
        for gene in kegg_cog:
            if (kegg_cog[gene]['KEGG'] != 'unknown'):
                kos = kegg_cog[gene]['KEGG'].split(';')
                if len(kos) == 1:
                    valid_genes.append(gene)
        valid_genes = list(set(valid_genes).intersection(set(df.columns)))
        df = df[valid_genes]
        print (str(datetime.now()) + " - EM shape: " + str(df.shape))
        return df
    
    

class EMGenesBuilder(DataBuilder):
    code_version = 1
    _composingitems = ()
    name = "EMGenesBuilderRaw2"
    sub_builders = SubDataBuilder('microbiome', \
                        (
                        ), \
                        (
                         EMGenesConcat(),
                         FixStuff(),
                         BuildKOdfTransformer(),
                        FDSelect(GenotekHandling.Without, StudyTypeHandling.PNP1AndBread, RemainingHandling.TakeLarger),
                        MinReadFilter('All', MetagenomeMDExposer.F_PostHGF, 1e7, remove = True),\
                         
                         RemoveGenesWithMoreThanOneKOTransfomer(),
                        RemoveLowAbundantGenesTransformer(),
                        UniteCorrelatedEMGenesWithinKOTransformer()
                         
                         )),
    _prefilters = ()

        


class MetabolonEMPairwiseCorrelationWriter(SingleExplorer):
    print (str(datetime.now()) + " - MetabolonEMPairwiseCorrelationWriter Class is parsed")
#     ppmet, samps, metabs = get_preprocessed()
    with open('/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Metabolon/tmp_files/ppmet.df', 'rb') as handle:
        ppmet = pickle.load(handle)
    with open('/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Metabolon/tmp_files/samps.df', 'rb') as handle:
        samps = pickle.load(handle)
    with open('/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Metabolon/tmp_files/metabs.df', 'rb') as handle:
        metabs = pickle.load(handle)
    
    
    b_A = AdHocBuilder(ppmet)
    print (str(datetime.now()) + " - EMGenesBuilder()")
    b_B = EMGenesBuilder()
    builder_A = b_A
    builder_B = b_B
    dropfeatures_A = ()
    dropfeatures_B = ()
#     corr_method = staticmethod(intersected_v_method_mincount(spearmanr, 200))
#     corr_method_wrapper = None
#     ignore_zero = False
#     fdr_per_column = False
#     outpath = '/net/mraid08/export/jafar/Microbiome/Analyses/Metabolon2/tmp'
#     corr_panel = {}
#     fdr_threshold = 0.05
    name = 'Itzak_KO'
    def internal_analysis(self, A, B):
        
        
        print (str(datetime.now()) + " - Start internal analysis")
#         B.T.to_csv(EM_GENES_PATH_READY)
        
#         if os.path.exists(EM_GENES_PATH_READY):
#         if EM_GENES_PATH_10 == EM_GENES_PATH_RAW:
#             print (str(datetime.now()) + " - Load csv")
#             self.B = read_csv(EM_GENES_PATH_READY)
#             print (str(datetime.now()) + " - Transposing csv")
#             self.B = self.B.set_index('Unnamed: 0').T #TODO: check if this helps
#             self.B = concat((self.B.reset_index(drop = True), \
#                         self.B.reset_index()['index'].str.split('---').apply(Series).\
#                             rename(columns = {0:'CLIENT_IDENTIFIER'})), axis = 1)
# #             self.B['ConnectionID'] = self.B.ConnectionID.astype(int)
#             self.B = self.B.set_index('CLIENT_IDENTIFIER')
#         else:
# #             assert 1 == 0
# #             print (str(datetime.now()) + " - Transposing csv")
        self.B = B
        
        self.A = A.reset_index(0, True).loc[self.B.index]
        self.A = self.A.iloc[:,[i for i in range(self.A.shape[1]) if len(self.A.iloc[:,i].unique()) > 2]]
        
#         print self.A.shape
#         with open('/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Metabolon/tmp_files/ppmt_483x1238.df', 'wb') as handle:
#             pickle.dump(self.A, handle, protocol=pickle.HIGHEST_PROTOCOL)
        print self.A.shape
        print self.B.shape
        print (str(datetime.now()) + " - End internal analysis")
        
M = MetabolonEMPairwiseCorrelationWriter()
M.run()