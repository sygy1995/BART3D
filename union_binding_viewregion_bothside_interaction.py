import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect

import utils

# === @marvinquiet remove normalization
def write_out_binding_interactions_all_chroms(binding_df,data,viewregion,resolution,flag,outdir):
    
    binding_data_out = open(outdir+os.sep+'{}_res{}_{}.csv'.format(data,resolution,flag),'w')
    columns = np.arange(-1*viewregion+resolution,viewregion,resolution)
    binding_data_out.write('{}\t{}\n'.format('id','\t'.join(map(str,columns))))
            
    for chrom in utils.chroms:
        matrix_file = '{}/{}_{}_res_{}_view_region_{}.csv'.format(outdir,data,chrom,resolution, viewregion)
        if os.path.isfile(matrix_file):
            with open(matrix_file) as matrix_inf:
                matrix_df = pd.read_csv(matrix_inf,sep='\t',index_col=0)
        
            # normalization of the data, each interaction score is divided by average interaction score with same distances
            matrix_df = matrix_df/matrix_df.mean()
            matrix_df = utils.mirror_concat(matrix_df)
            binding_df_chr = binding_df.loc[binding_df['chr']==chrom]
            for binding_id in binding_df_chr.index:
                binding_pos = binding_df.loc[binding_id].middle
                view_pos = binding_pos//resolution*resolution
                view_pos_binding_data = matrix_df.loc[view_pos][map(str,columns)].values  
                view_pos_binding_data = np.round(view_pos_binding_data,2)
                binding_data_out.write('{}\t{}\n'.format(binding_id,'\t'.join(map(str,view_pos_binding_data))))
    binding_data_out.close()
 
        
def main(view_region,normalization,chrom):
    
    outdir = 'f2_union_binding_view2M_bothside_interaction'
    os.makedirs(outdir,exist_ok=True)
    
    viewregion = 2000000
    
    for normalization in ['raw','iced']:
        for data in ['Jurkat','A6010','Cutll1','PD31','PD9']:
            for resolution in [5000,10000,50000,100000]:          
                write_out_binding_interactions_all_chroms(all_bindings,data,normalization,viewregion,resolution,'union',outdir)







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
