import os, sys
import argparse

from read_normalize_interaction import get_normalized_viewpoint_interaction
from get_interaction_change import compare_hic_interaction
import pandas as pd
import utils

import time

def main(args):
    os.makedirs(args.outdir, exist_ok=True)

    tick = time.time()
    # get resolution of Hi-C data
    res_c = utils.get_resolution(args.c_index)
    res_t = utils.get_resolution(args.t_index)
    if res_c != res_t :
        sys.stderr.write("Error: Resolutions from two Hi-C index file are not the same!\n")
        sys.exit(1)
    args.resolution = res_c

    # TODO: according to different species, get chrom and length for the species
    # args.species
    if args.species=='hg38':
        chroms = utils.chroms_hg38
    elif args.species=='mm10':
        chroms = utils.chroms_mm10
    else:
        sys.stderr.write("Error: Species has to be either hg38 or mm10.\n")
        sys.exit(1)

    # get prefix of file
    c_prefix = os.path.splitext(os.path.basename(args.c_index))[0]
    t_prefix = os.path.splitext(os.path.basename(args.t_index))[0]

    # get matrix dataframe
    matrix_df_control = pd.read_csv(args.c_matrix,sep='\t',header=None)
    matrix_df_control.columns = ['id1','id2','count']
    matrix_df_treatment = pd.read_csv(args.t_matrix,sep='\t',header=None)
    matrix_df_treatment.columns = ['id1','id2','count']

    # get index dict
    index = utils.get_index(args.c_index)

    # initiaze output data
    output_file_name = args.outdir+os.sep+'{}_over_{}_res{}_view{}.bed'.format(t_prefix,c_prefix,args.resolution,args.region)
    compr_data_out = open(output_file_name,'w')

    # the rest of operations are divided to each chromosomes
    for chrom in chroms:
        sys.stdout.write("Getting normalized Hi-C viewpoint interaction data from Juicer format for chromosome {}..\n".format(chrom))

        interaction_normalized_control_np = get_normalized_viewpoint_interaction(chrom,index,matrix_df_control,args.region,args.resolution,args.species)
        interaction_normalized_treatment_np = get_normalized_viewpoint_interaction(chrom,index,matrix_df_treatment,args.region,args.resolution,args.species)

        sys.stdout.write("Getting Hi-C viewpoint interaction change for chromosome {}..\n".format(chrom))

        compare_hic_interaction(interaction_normalized_control_np,interaction_normalized_treatment_np,args.resolution,output_file_name,chrom,args.species)

    tock = time.time()-tick
    print('finished '+str(tock))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Identify transcription factors which may cause the different interactions between two Hi-C profiles.')

    parser.add_argument('-ci', '--controlindex', action='store', type=str, dest='c_index', required=True, help='Hi-C index file for control', metavar='<str>')
    parser.add_argument('-cm', '--controlmatrix', action='store', type=str, dest='c_matrix', required=True, help='Hi-C matrix file for control', metavar='<str>')
    parser.add_argument('-ti', '--treatindex', action='store', type=str, dest='t_index', required=True, help='Hi-C index file for treat', metavar='<str>')
    parser.add_argument('-tm', '--treatmatrix', action='store', type=str, dest='t_matrix', required=True, help='Hi-C matrix file for treat', metavar='<str>')

    parser.add_argument('-s', '--species',dest='species',type=str,metavar='<species>',choices=['hg38','mm10'],required=True,help='Species, please choose from "hg38" or "mm10".')

    parser.add_argument('-o', '--outdir', action='store', type=str, dest='outdir', help='output directory for Hi-C bart', metavar='<str>', default='hic_interaction_output/')
    parser.add_argument('-r', '--region', action='store', type=int, dest='region', help='Regions to expand when finding interactions', metavar='<int>', default=500000)

    args = parser.parse_args()
    if(len(sys.argv)) < 5: # two index, two matrix file, one species, are required
        parser.print_help()
        sys.exit(1)

    main(args)
