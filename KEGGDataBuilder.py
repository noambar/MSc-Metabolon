###################################################################################################
# File: KEGGDataBuilder.py
# Version: 0.0
# Date: 14.8.2017
# Noam Bar, noam.bar@weizmann.ac.il
#
# 
# Python version: 2.7
###################################################################################################

import Utils
import urllib2
from datetime import datetime
import sys
import re
import os
from queue.qp import qp
from addloglevels import sethandlers

KEGG_BASE = 'http://rest.kegg.jp/'
KEGG_COG_DICT_PATH = '/net/mraid08/export/jafar/Microbiome/Data/Databases/GeneSet2014/Annotation/KEGG_COG_dict.dat'
KEGG_DB_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/KEGG_DB_12_12_2017'
base_urls = {'protein':'http://www.genome.jp/dbget-bin/www_bget?-f+-n+a+', 
                 'dna':'http://www.genome.jp/dbget-bin/www_bget?-f+-n+n+'}
KEGG_ORGANISM_LIST_URL = 'http://rest.kegg.jp/list/organism'
url_organism_list = urllib2.urlopen(KEGG_ORGANISM_LIST_URL, timeout=20).readlines()
url_organism_dict = {org.split()[1]:org.split('\n')[0] for org in url_organism_list}
# This class holds KEGG links between different components.
# It uses KEGG public data from http://rest.kegg.jp/
class KEGG_data_holder(object):
    
    def __init__(self, keep_unique_bacgenes = True):
        """Constructor"""
        self._build_urls()    
        self._build_dicts()
        self._build_KEGG_COG_dict(keep_unique_bacgenes)
        self._link_dicts('compound_pathway', 'pathway_ko')
        self._link_dicts('compound_reaction', 'reaction_ko')
        self._link_dicts('compound_module', 'module_ko')
        self._link_dicts('drug_pathway', 'pathway_ko')
        self._reverse_dicts('compound_pathway_ko')
        self._reverse_dicts('compound_module_ko')
        self._reverse_dicts('compound_reaction_ko')
        
           
    def _build_urls(self):
        """Define link's url."""
        links = [('compound', 'pathway'), ('compound', 'reaction'), ('compound', 'module'),
                 ('drug', 'pathway'), ('drug', 'ko'), ('reaction', 'ko'), ('pathway', 'ko'),
                 ('module', 'ko')]
        self.urls = {'_'.join([link[0], link[1]]): KEGG_BASE + 
                     '/'.join(['link', link[1], link[0]]) for link in links}

    
    def _build_KEGG_COG_dict(self, keep_unique_bacgenes, KEGG_COG_dict_path = KEGG_COG_DICT_PATH):
        """Build a dictionary which relates KOs to bacterial genes."""
        
        print(str(datetime.now()) + 
              " - [KEGG_data_holder._build_KEGG_COG_dict] loading KEGG_COG dict")
        self.dicts['KEGG_COG'] = Utils.Load(KEGG_COG_dict_path)
        self.dicts['ko_bacgene'] = {}
        print(str(datetime.now()) + 
              " - [KEGG_data_holder._build_KEGG_COG_dict] parsing ko_bacgene dict")
        for gene in self.dicts['KEGG_COG']:
            if (self.dicts['KEGG_COG'][gene]['KEGG'] != 'unknown'):
                kos = self.dicts['KEGG_COG'][gene]['KEGG'].split(';')
                if keep_unique_bacgenes:
                    if len(kos) > 1:
                        continue
                for ko in kos:
                    kor = 'ko:' + ko
                    if (kor in self.dicts['ko_bacgene']):
                        self.dicts['ko_bacgene'][kor] += [gene]
                    else:
                        self.dicts['ko_bacgene'][kor] = [gene]
        self.dicts['cog_bacgene'] = {}
        print(str(datetime.now()) + 
              " - [KEGG_data_holder._build_KEGG_COG_dict] parsing cog_bacgene dict")
        for gene in self.dicts['KEGG_COG']:
            if (self.dicts['KEGG_COG'][gene]['COG'] != 'unknown'):
                cogs = self.dicts['KEGG_COG'][gene]['COG'].split(';')
                if keep_unique_bacgenes:
                    if len(cogs) > 1:
                        continue
                for cog in cogs:
                    cogr = 'cog:' + cog
                    if (cogr in self.dicts['cog_bacgene']):
                        self.dicts['cog_bacgene'][cogr] += [gene]
                    else:
                        self.dicts['cog_bacgene'][cogr] = [gene]
      
    def _build_dicts(self):
        """Build all dictionaries"""
        self.dicts = {}
        for url in self.urls:
            print(str(datetime.now()) + 
                  " - [KEGG_data_holder._build_dicts] building url " + url)
            n_tries = 0
            while (n_tries < 10):
                try:
                    self.dicts[url] = self._KEGG_link_to_dict(urllib2.urlopen(self.urls[url], 
                                                                              timeout=20).readlines())
                    break
                except:
                    n_tries += 1
            if (n_tries == 10):
                sys.stderr.write("Couldn't build url " + url + "\n")
        
        
    def _KEGG_link_to_dict(self, KEGG_link):
        """Convert KEGG link format to dictionary."""
        KEGG_dict = {}
        for line in KEGG_link:
            k_v = line.split('\n')[0].split('\t')
            if (k_v[0] not in KEGG_dict):
                KEGG_dict[k_v[0]] = [k_v[1]]
            else:
                KEGG_dict[k_v[0]] += [k_v[1]]
        return KEGG_dict
    
    def _link_dicts(self, dict1, dict2):
        """Given two dictionaries entries, link first element of first dict to last element 
        of second dict"""
        print(str(datetime.now()) + " - [KEGG_data_holder._link_dicts] linking " + 
              dict1 + ", " + dict2)
        if (dict1.split('_')[-1] != dict2.split('_')[0]):
            sys.stderr.write("ERROR [KEGG_data_holder._link_dicts]: Please enter two dictionaries \
            such that the values of the first are the keys of the second.\n")
            return
        new_dict_name = dict1 + '_' + '_'.join(dict2.split('_')[1:])
        self.dicts[new_dict_name] = {}
        for k in self.dicts[dict1]:
            final_vals = []
            for v in self.dicts[dict1][k]:
                if (v in self.dicts[dict2]):
                    final_vals += self.dicts[dict2][v]
            self.dicts[new_dict_name][k] = final_vals
            
    def _reverse_dicts(self, dict_to_reverse):
        if dict_to_reverse not in self.dicts:
            return
        names = dict_to_reverse.split('_')
        names.reverse()
        new_name = '_'.join(names)
        self.dicts[new_name] = {}
        for k1 in self.dicts[dict_to_reverse]:
            for v in self.dicts[dict_to_reverse][k1]:
                if v in self.dicts[new_name]:
                    self.dicts[new_name][v].append(k1)
                else:
                    self.dicts[new_name][v] = [k1]
        for k2 in self.dicts[new_name]:
            self.dicts[new_name][k2] = list(set(self.dicts[new_name][k2]))
        return
    
    def _merge_dicts(self, prefix_dict, suffix_dict):
        if '_'.join([prefix_dict, suffix_dict]) in self.dicts:
            return
        new_dict_name = '_'.join([prefix_dict, suffix_dict])
        self.dicts[new_dict_name] = {}
        for other_dict in self.dicts.keys():
            dict_elements = other_dict.split('_')
            if dict_elements[0] == prefix_dict and dict_elements[-1] == suffix_dict:
                for k in self.dicts[other_dict]:
                    if k in self.dicts[new_dict_name]:
                        self.dicts[new_dict_name][k] += self.dicts[other_dict][k]
                    else:
                        self.dicts[new_dict_name][k] = self.dicts[other_dict][k]
        for k in self.dicts[new_dict_name]:
            self.dicts[new_dict_name][k] = list(set(self.dicts[new_dict_name][k]))
        return
    
    def get_dicts(self):
        return self.dicts
    
    def get_link_dicts(self, dict1, dict2, key1):
        """Given two dictionaries and a key for the first one, return all values chained"""
        final_vals = []
        if (key1 not in self.dicts[dict1]):
#             sys.stderr.write("ERROR [KEGG_data_holder.get_link_dicts]: Key not in " + dict1 + ".\n")
            return final_vals
        possible_keys = self.dicts[dict1][key1]
        for k in possible_keys:
            if (k in self.dicts[dict2]):
                final_vals += self.dicts[dict2][k]
        return final_vals
    
def KEGG_read_all_KOs_from_web():
    KEGG_KO_LIST_URL = 'http://rest.kegg.jp/list/ko'
    url_kos_list = urllib2.urlopen(KEGG_KO_LIST_URL, timeout=20).readlines()
    url_kos_dict = {ko.split('\t')[0]:ko.split('\n')[0].split('\t')[1] for ko in url_kos_list}
    Utils.Write(KEGG_DB_DIR + '/kos_dict.dat', url_kos_dict)
    return

    
def KEGG_create_gene_db_from_KO(ko, outdir, aa_dna='protein'):
    base_urls = {'protein':'http://www.genome.jp/dbget-bin/www_bget?-f+-n+a+', 
                 'dna':'http://www.genome.jp/dbget-bin/www_bget?-f+-n+n+'}
    KEGG_ORGANISM_LIST_URL = 'http://rest.kegg.jp/list/organism'
    url_organism_list = urllib2.urlopen(KEGG_ORGANISM_LIST_URL, timeout=20).readlines()
    url_organism_dict = {org.split()[1]:org.split('\n')[0] for org in url_organism_list}
    LEGAL_SEQ = ['protein', 'dna']
    if aa_dna not in LEGAL_SEQ:
        print(str(datetime.now()) + " - [KEGGDataBuilder.KEGG_create_gene_db_from_KO] aa_dna must be one of: " 
              + ", ".join(LEGAL_SEQ))
        return
    BASE_URL = base_urls[aa_dna]
    
    if not re.match('^K[0-9]+$', ko):
        print(str(datetime.now()) + " - [KEGGDataBuilder.KEGG_create_gene_db_from_KO] KO name " 
              + ko + " isn't legal")
        return
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    url = 'http://rest.kegg.jp/link/genes/' + ko
    url_genes = urllib2.urlopen(url, timeout=20).readlines()
    genes = [gene.split('\n')[0].split()[1] for gene in url_genes]
    genes_seq = {}
    with open('/'.join([outdir, ko + '_'+ aa_dna + '.fa']), 'w') as handle:
        print(str(datetime.now()) + " - [KEGGDataBuilder.KEGG_create_gene_db_from_KO] will parse " 
              + str(len(genes)) + " genes")
        c = 0
        for gene in genes:
            c += 1
            if c % 10 == 0:
                print c
            if not bool(re.match('.*Bacteria.*', url_organism_dict[gene.split(':')[0]])):
                print gene.split(':')[0]
                continue
            url = BASE_URL + gene
            url_text = urllib2.urlopen(url, timeout=20).readlines()
            text = ''.join([t.split('\n')[0] for t in url_text if (t[0] != '<') & bool((re.match('[a-zA-Z]+', t)))])
            if text in genes_seq.values():
                continue
            genes_seq[gene] = text
            handle.write('>' + gene + '\n')
            handle.write(text + '\n')
            
    return  
    

    
        

# kegg = KEGG_data_holder()
# kegg._link_dicts('compound_pathway', 'pathway_ko')
# kegg._link_dicts('compound_reaction', 'reaction_ko')
# kegg._link_dicts('compound_module', 'module_ko')
# kegg._link_dicts('drug_pathway', 'pathway_ko')
# kegg._link_dicts('compound_ko', 'ko_bacgene') # too long to run