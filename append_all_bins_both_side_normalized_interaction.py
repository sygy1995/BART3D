import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect

import utils


# === @marvinquiet: remove normalization
def write_out_binding_interactions_sep_chroms(data,viewregion,resolution,outdir,species):
    columns = np.arange(-1*viewregion+resolution,viewregion,resolution)
    
    if species=='hg38':
        chroms = utils.chroms_hg38
    elif species=='mm10':
        chroms = utils.chroms_mm10

    for chrom in chroms:
        matrix_file = outdir+os.sep+'{}_{}_res_{}_view_region_{}.csv'.format(data,chrom,resolution, viewregion)
        all_region_data_out = outdir+os.sep+'{}_res{}_view{}_{}.csv'.format(data,resolution,viewregion,chrom)

        if os.path.isfile(matrix_file):
            with open(matrix_file) as matrix_inf:
                matrix_df = pd.read_csv(matrix_inf,sep='\t',index_col=0)
                #####
                # process the matrix df
                #####
        
            # normalization of the data, each interaction score is divided by average interaction score with same distances
            matrix_df = matrix_df/matrix_df.mean()
            matrix_df = utils.mirror_concat(matrix_df)

            # === @marvinquiet: change the column index data type to int
            matrix_df.columns = matrix_df.columns.astype(int)
            matrix_df = matrix_df[[i for i in map(int,columns)]]
            # === end

            matrix_df.fillna(0)
            matrix_df.to_csv(all_region_data_out)
    
        
        
def main(view_region,normalization,chrom):
    
    outdir = 'f1_append_2M_BS_normalized_interaction_sep_chr'
    os.makedirs(outdir,exist_ok=True)
    
    viewregion = 1000000
    
    for normalization in ['raw','iced']:
        for data in ['Jurkat','A6010','Cutll1','PD31','PD9']:
            #for resolution in [5000,10000,50000,100000]:    
            for resolution in [5000]:          
                write_out_binding_interactions_sep_chroms(data,normalization,viewregion,resolution,outdir)




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--viewregion', action = 'store', type = int,dest = 'viewregion', help = 'input file of', metavar = '<int>')
    parser.add_argument('-n', '--normalization', action = 'store', type = str,dest = 'normalization', help = 'input file of', metavar = '<str>')
    parser.add_argument('-c', '--chrom', action = 'store', type = str,dest = 'chrom', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.viewregion,args.normalization,args.chrom)
