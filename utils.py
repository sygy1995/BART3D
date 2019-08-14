import os
import numpy as np
import pandas as pd

##########################################
## basic models for commonly used data 
##########################################
current_dir = os.path.dirname(os.path.realpath(__file__))
chrom_size_file_hg38 = os.path.join(current_dir,'hg38_clean.chrom.sizes')
chrom_size_df_hg38 = pd.read_csv(chrom_size_file_hg38,sep='\t',header=None,index_col=0);
chrom_size_df_hg38.columns = ['len']

chroms_hg38 = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']

chrom_size_file_mm10 = os.path.join(current_dir,'mm10_clean.chrom.sizes')
chrom_size_df_mm10 = pd.read_csv(chrom_size_file_mm10,sep='\t',header=None,index_col=0);
chrom_size_df_mm10.columns = ['len']

chroms_mm10 = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
             'chr18','chr19','chrX','chrY']


##########################################
## basic functions for build mirror for matrix 
##########################################

# get resolution from index file: chr, start, end, index
def get_resolution(hic_index_file):
    findex = open(hic_index_file, 'r')
    line = findex.readline().strip().split()
    resolution = int(line[2])-int(line[1])
    findex.close()
    return resolution

def get_index(order_index_file):
    f = open(order_index_file)
    index = {}
    for line in f:
        sline = line.strip().split()
        if sline[0] not in index:
            index[sline[0]] = {}
        index[sline[0]][int(sline[3])] = int(sline[1])
    return index


