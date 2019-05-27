import os
import numpy as np
import pandas as pd

chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',\
             'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',\
             'chr18','chr19','chr20','chr21','chr22','chrX','chrY']

##########################################
## basic models for commonly used data 
##########################################
chrom_size_file = 'hg38_clean.chrom.sizes'
chrom_size_df = pd.read_csv(chrom_size_file,sep='\t',header=None,index_col=0);
chrom_size_df.columns = ['len']

# get all bindings, TODO: canbe written as a file
union_CTCF_kept_ids_sample_thre_2_file = 'union_CTCF_kept_ids_sample_thre_2.csv'
all_binding_middle_domain_file = 'all_CTCF_domainInfo.csv'
with open(union_CTCF_kept_ids_sample_thre_2_file) as filter_inf, open(all_binding_middle_domain_file) as all_inf:
    kept_ids = [int(i.strip()) for i in filter_inf.readlines()]
    all_bindings = pd.read_csv(all_inf,sep='\t')
all_bindings = all_bindings.loc[kept_ids]


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

# get matrix df according to view region, default=500,000bps
def get_matrix_df(matrix_file, view_region):
    if os.path.getsize(matrix_file) == 0:
        return None

    with open(matrix_file) as matrix_inf:
        matrix_df = pd.read_csv(matrix_inf,sep='\t',header=None)
        matrix_df.columns = ['x','y','score']
    matrix_df = matrix_df[matrix_df['y']-matrix_df['x']<view_region]
    # get the distances between x and y
    matrix_df['dis'] = matrix_df['y']-matrix_df['x']
    return matrix_df

