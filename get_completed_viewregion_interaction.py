import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect

import utils

# @marvinquiet: remove normalization
# def write_out_interactions(matrix_df,view_region,outdir,data,resolution,chrom):
#     # for each position/view_pos in chrom, get the interaction with downstream loci within view_region
#     outfile_name = outdir+os.sep+'{}_{}_res_{}_view_region_{}.csv'.format(data,chrom,resolution,view_region)
#     if os.path.isfile(outfile_name):
#         print('Do not re-write the files')
#         exit()
#     out_file = open(outfile_name,'w')
#     view_poses = np.arange(0,utils.chrom_size_df.loc[chrom,'len'],resolution)

#     out_file.write('{}\t{}\n'.format('dis','\t'.join(map(str,np.arange(0,view_region,resolution)))))
#     for view_pos in view_poses:
#         view_df = matrix_df[matrix_df['x']==view_pos]
#         view_df = view_df[['score','dis']]; 
#         standard_df = pd.DataFrame(columns = ['dis'],index = np.arange(0,view_region,resolution))
#         standard_df['dis']=standard_df.index
#         # the matrix_df has been restricted within view-region, here we could use outer or left
#         # here is to get each of the interaction scores within the view region, if the interaction score exist in view_df, fill with 0
#         merge_df = pd.merge(standard_df,view_df,how='outer').fillna(0)
#         out_file.write('{}\t{}\n'.format(view_pos,'\t'.join(map(str,merge_df.score.values))));
#     out_file.close()

def write_out_interactions(matrix_dict,view_region,outdir,data,resolution,chrom):
    # for each position/view_pos in chrom, get the interaction with downstream loci within view_region
    outfile_name = outdir+os.sep+'{}_{}_res_{}_view_region_{}.csv'.format(data,chrom,resolution,view_region)
    if os.path.isfile(outfile_name):
        print('Do not re-write the files')
        exit()
    out_file = open(outfile_name,'w')
    view_poses = np.arange(0,utils.chrom_size_df.loc[chrom,'len'],resolution)

    out_file.write('{}\t{}\n'.format('dis','\t'.join(map(str,np.arange(0,view_region,resolution)))))
    for view_pos in view_poses:
        if view_pos in matrix_dict:
            interactions = [0.0]*int(view_region/resolution)
            for anchor in matrix_dict[view_pos]:
                interactions[int((anchor-view_pos)/resolution)] = matrix_dict[view_pos][anchor]
        else:
            interactions = [0.0]*int(view_region/resolution)
        interactions.insert(0,view_pos.item())
        out_file.write('\t'.join(map(str,interactions))+'\n')
    


def main(chrom):
    
    indir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/11_HiC_analysis/f1_preprocess/union_binding_processed_data/panos/transformed_matrix'
    outdir = 'f1_viewpoint_2M_interaction_abstraction'
    os.makedirs(outdir,exist_ok=True)
    view_region = 2000000

    for data in ['Jurkat','A6010','Cutll1','PD31','PD9']:
        for resolution in [5000,10000,50000,100000]:
            for normalization in ['raw','iced']:
                matrix_file = indir+os.sep+'{}/{}/{}_{}_{}_{}.matrix'.format(data,resolution,data,resolution,normalization,chrom)
                with open(matrix_file) as matrix_inf:
                    matrix_df = pd.read_csv(matrix_inf,sep='\t',header=None)
                    matrix_df.columns = ['x','y','score'.format(data)]
                matrix_df = matrix_df[matrix_df['y']-matrix_df['x']<view_region]
                # get the distances between x and y
                matrix_df['dis'] = matrix_df['y']-matrix_df['x']
                #print(matrix_df)
                write_out_interactions(matrix_df,view_region,outdir,data,resolution,normalization,chrom)
                #exit()
    




if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    #parser.add_argument('-v', '--viewregion', action = 'store', type = int,dest = 'viewregion', help = 'input file of', metavar = '<int>')
    #parser.add_argument('-d', '--data', action = 'store', type = str,dest = 'data', help = 'input file of', metavar = '<int>')
    #parser.add_argument('-n', '--normalization', action = 'store', type = str,dest = 'normalization', help = 'input file of', metavar = '<str>')
    parser.add_argument('-c', '--chrom', action = 'store', type = str,dest = 'chrom', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.chrom)
