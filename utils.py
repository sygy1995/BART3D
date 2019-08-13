import os
import numpy as np
import pandas as pd

##########################################
## basic models for commonly used data 
##########################################
chrom_size_file_hg38 = '/home/yifan/Documents/hic_bart/hg38_clean.chrom.sizes'
chrom_size_df_hg38 = pd.read_csv(chrom_size_file_hg38,sep='\t',header=None,index_col=0);
chrom_size_df_hg38.columns = ['len']

chroms_hg38 = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']

chrom_size_file_mm10 = '/home/yifan/Documents/hic_bart/mm10_clean.chrom.sizes'
chrom_size_df_mm10 = pd.read_csv(chrom_size_file_mm10,sep='\t',header=None,index_col=0);
chrom_size_df_mm10.columns = ['len']

chroms_mm10 = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
             'chr18','chr19','chrX','chrY']


##########################################
## basic functions for build mirror for matrix 
##########################################
# mirror the matrix
def mirror_concat(df):
    for column_pos in np.arange(1,len(df.columns)):
        column = df.columns[column_pos]
        mir_column = '-{}'.format(column)
        new_values = np.append([0]*column_pos,df[column].values)[:len(df.index)]
        df[mir_column]=new_values
    return df

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

# get matrix df according to view region, default=500,000bps
# def get_matrix_df(matrix_file, view_region):
#     if os.path.getsize(matrix_file) == 0:
#         return None

#     with open(matrix_file) as matrix_inf:
#         matrix_df = pd.read_csv(matrix_inf,sep='\t',header=None)
#         matrix_df.columns = ['x','y','score']
#     matrix_df = matrix_df[matrix_df['y']-matrix_df['x']<view_region]
#     # get the distances between x and y
#     matrix_df['dis'] = matrix_df['y']-matrix_df['x']
#     return matrix_df

def get_matrix_dict(matrix_file, view_region):
    if os.path.getsize(matrix_file) == 0:
        return None
    f = open(matrix_file)
    matrix_dict = {}
    for line in f:
        sline = line.strip().split()
        start = int(sline[0])
        end = int(sline[1])
        count = float(sline[2])
        if (end-start) < view_region:
            if start not in matrix_dict:
                matrix_dict[start] = {}
            matrix_dict[start][end] = count
    return matrix_dict



