###################################################################################################
# File: RunDiamond.py
# Version: 0.0
# Date: 6.12.2017
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
import time
import urllib2, sys
from time import sleep
# from BeautifulSoup import BeautifulSoup

KEGG_BASE = 'http://rest.kegg.jp/'
KEGG_COG_DICT_PATH = '/net/mraid08/export/jafar/Microbiome/Data/Databases/GeneSet2014/Annotation/KEGG_COG_dict.dat'
KEGG_DB_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/KEGG_DB_12_12_2017'
# base_urls = {'protein':'http://www.genome.jp/dbget-bin/www_bget?-f+-n+a+', 
#                  'dna':'http://www.genome.jp/dbget-bin/www_bget?-f+-n+n+'}
base_urls = {'protein':'/aaseq', 'dna':'/ntseq'}
KEGG_ORGANISM_LIST_URL = 'http://rest.kegg.jp/list/organism'
url_organism_list = urllib2.urlopen(KEGG_ORGANISM_LIST_URL, timeout=20).readlines()
url_organism_dict = {org.split()[1]:org.split('\n')[0] for org in url_organism_list}

def KEGG_read_all_KOs_from_web():
    KEGG_KO_LIST_URL = 'http://rest.kegg.jp/list/ko'
    url_kos_list = urllib2.urlopen(KEGG_KO_LIST_URL, timeout=20).readlines()
    url_kos_dict = {ko.split('\t')[0]:ko.split('\n')[0].split('\t')[1] for ko in url_kos_list}
    Utils.Write(KEGG_DB_DIR + '/kos_dict.dat', url_kos_dict)
    return

def KEGG_build_KO_genes_dict():
    if not os.path.exists(KEGG_DB_DIR + '/kos_dict.dat'):
        KEGG_read_all_KOs_from_web()
    kos_dict = Utils.Load(KEGG_DB_DIR + '/kos_dict.dat')
    kos = kos_dict.keys()
    kos_genes_dict = {}
    c=0
    for ko in kos:
        c+=1
        if c % 100 == 0:
            print datetime.now()
            Utils.Write(KEGG_DB_DIR + '/kos_genes_dict.dat', kos_genes_dict)
            print c
        url = 'http://rest.kegg.jp/link/genes/' + ko.split('ko:')[1]
        url_genes = urllib2.urlopen(url, timeout=20).readlines()
        try:
            genes = [gene.split('\n')[0].split()[1] for gene in url_genes]
        except:
            print ko
        kos_genes_dict[ko] = genes
    Utils.Write(KEGG_DB_DIR + '/kos_genes_dict.dat', kos_genes_dict)
    return

def KEGG_build_genes_KO_dict():
    if not os.path.exists(KEGG_DB_DIR + '/kos_dict.dat'):
        KEGG_read_all_KOs_from_web()
    if not os.path.exists(KEGG_DB_DIR + '/kos_genes_dict.dat'):
        print "KO - genes file doesn't exist. Run: KEGG_build_KO_genes_dict()"
        return
    kos_genes_dict = Utils.Load(KEGG_DB_DIR + '/kos_genes_dict.dat')
    kos = kos_genes_dict.keys()
    genes_kos_dict = {}
    c=0
    for ko in kos:
        c+=1
        if c % 100 == 0:
            print datetime.now()
            print c
        genes = kos_genes_dict[ko]
        for g in genes:
            if g in genes_kos_dict:
                genes_kos_dict[g].append(ko)
            else:
                genes_kos_dict[g] = [ko]

    Utils.Write(KEGG_DB_DIR + '/genes_kos_dict.dat', genes_kos_dict)
    return
     
def KEGG_download_and_create_fasta_file_from_genes(genes_path, taxa = 'Bacteria', aa_dna = 'protein', genes_per_batch = 10000):
    KEGG_TEMP_FASTA = KEGG_DB_DIR + '/temp_fasta_' + aa_dna + '/'
    BASE_URL = base_urls[aa_dna]
    if not os.path.exists(KEGG_TEMP_FASTA): os.makedirs(KEGG_TEMP_FASTA)
          
    if not os.path.exists(KEGG_DB_DIR + '/kos_genes_dict.dat'):
        print "KO - genes file doesn't exist. Run: KEGG_build_KO_genes_dict()"
        return
    genes = [gene.split('\n')[0] for gene in open(genes_path, 'r').readlines()]
#     kos_genes_dict = Utils.Load(KEGG_DB_DIR + '/kos_genes_dict.dat')
    with qp(jobname = 'KGG2fasta', q=['himem7.q'], mem_def = '1G', trds_def = 1, tryrerun=True, 
            max_u = 210, delay_batch=15) as q:
        os.chdir("/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/temp_q_dir/")
        q.startpermanentrun()
        upload_create_fasta_files_jobs(q, KEGG_TEMP_FASTA, genes, BASE_URL, taxa, genes_per_batch)

    output_fasta = KEGG_DB_DIR + '_'.join(['/KEGG_genes', taxa, aa_dna]) + '.fa'
    temp_files = os.listdir(KEGG_TEMP_FASTA)
    temp_files.sort()
    
    with open(output_fasta, 'w') as handle:
        for temp_file in temp_files:
            print temp_file
            fasta_content = open(KEGG_TEMP_FASTA + '/' + temp_file, 'r').readlines()
            for line in fasta_content:
                handle.write(line)
            os.remove(KEGG_TEMP_FASTA + '/' + temp_file)
    os.removedirs(KEGG_TEMP_FASTA)
    return

def upload_create_fasta_files_jobs(q, KEGG_TEMP_FASTA, genes, BASE_URL, taxa, genes_per_batch = 10000):
    waiton = []
    for i in range(0, len(genes), genes_per_batch):
        if i+genes_per_batch > len(genes):
            waiton.append(q.method(create_fasta_files, (KEGG_TEMP_FASTA, genes[i:], BASE_URL, taxa, i, genes_per_batch)))
        else:
            waiton.append(q.method(create_fasta_files, (KEGG_TEMP_FASTA, genes[i:i+genes_per_batch], BASE_URL, taxa, i, genes_per_batch)))
        
    res = q.waitforresults(waiton)
    return

def create_fasta_files(KEGG_TEMP_FASTA, genes, BASE_URL, taxa, start, genes_per_batch = 10000):
    with open('/'.join([KEGG_TEMP_FASTA, 'kegg_genes_' + str(start) + '.fa']), 'w') as handle:
        print start
        for i in range(0, len(genes), 5):
            end = i + 5
            if end > len(genes):
                end = len(genes)
            temp_genes = []
            for g in genes[i:end]:
                if g.split(':')[0] in url_organism_dict:
                    if bool(re.match('.*' + taxa + '.*', url_organism_dict[g.split(':')[0]])): temp_genes.append(g)
            if len(temp_genes) == 0:
                continue
#             print i
            time.sleep(0.5)
            url = 'http://rest.kegg.jp/get/' + '+'.join(temp_genes) + BASE_URL
            hdr = {'User-Agent': 'Mozilla/' + str(i)}
            try:
                req = urllib2.Request(url, headers=hdr)
                url_text = urllib2.urlopen(req).readlines()
                for l in url_text:
                    handle.write(l)
            except urllib2.HTTPError, e:
                print e.code
                print e.msg
                print i
#                 return
                os.remove('/'.join([KEGG_TEMP_FASTA, 'kegg_genes_' + str(start) + '.fa']))
                raise
            
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




def main():
    KEGG_download_and_create_fasta_file_from_genes('/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/KEGG_DB_12_12_2017/KEGG_allgenes.txt', 
                                                    taxa = 'Bacteria', aa_dna = 'protein', genes_per_batch = 10000)
    return

# if __name__ == "__main__":
#     sethandlers()
#     main()