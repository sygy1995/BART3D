import os, sys
import argparse

from BART3D.read_normalize_interaction import get_normalized_viewpoint_interaction
from BART3D.get_interaction_change import compare_hic_interaction
import pandas as pd
from BART3D import utils, score_on_UDHS, AUCcalc, StatTest, OptValidator



import time

def BART3D(options):
    # get library dirs
    args = OptValidator.opt_validate(options)


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
    output_file_name_up = args.outdir+os.sep+'{}_over_{}_res{}_view{}_upinteractions.bed'.format(t_prefix,c_prefix,args.resolution,args.region)
    output_file_name_down = args.outdir+os.sep+'{}_over_{}_res{}_view{}_downinteractions.bed'.format(t_prefix,c_prefix,args.resolution,args.region)
    compr_data_out_up = open(output_file_name_up,'w')
    compr_data_out_down = open(output_file_name_down,'w')

    # the rest of operations are divided to each chromosomes
    for chrom in chroms:
        sys.stdout.write("Getting normalized Hi-C viewpoint interaction data from HiC-Pro format for chromosome {}..\n".format(chrom))

        interaction_normalized_control_np = get_normalized_viewpoint_interaction(chrom,index,matrix_df_control,args.region,args.resolution,args.species)
        interaction_normalized_treatment_np = get_normalized_viewpoint_interaction(chrom,index,matrix_df_treatment,args.region,args.resolution,args.species)

        sys.stdout.write("Getting Hi-C viewpoint interaction change for chromosome {}..\n".format(chrom))

        compare_hic_interaction(interaction_normalized_control_np,interaction_normalized_treatment_np,args.resolution,output_file_name_up,output_file_name_down,chrom,args.species)

    # Region BART part
    sys.stdout.write('Start mapping the differential interaction profiles onto UDHS...\n')
    sys.stdout.flush()
    counting_up = score_on_UDHS.score_on_DHS(args,output_file_name_up)
    positions_up = sorted(counting_up.keys(),key=counting_up.get,reverse=True)
    counting_down = score_on_UDHS.score_on_DHS(args,output_file_name_down)
    positions_down = sorted(counting_down.keys(),key=counting_down.get,reverse=True)

    sys.stdout.write('BART Prediction starts...\n\nRank all DHS...\n')
    sys.stdout.flush()

    ofile_up = args.outdir+os.sep+'{}_over_{}_TF_upinteraction'.format(t_prefix,c_prefix)
    tf_aucs_up, tf_index_up = AUCcalc.cal_auc(args, positions_up,ofile_up)
    ofile_down = args.outdir+os.sep+'{}_over_{}_TF_downinteraction'.format(t_prefix,c_prefix)
    tf_aucs_down, tf_index_down = AUCcalc.cal_auc(args, positions_down,ofile_down)

    stat_file_up = ofile_up + '_bart_results.txt'
    stat_file_down = ofile_down + '_bart_results.txt'
    StatTest.stat_test(tf_aucs_up, tf_index_up, stat_file_up, args.normfile)
    StatTest.stat_test(tf_aucs_down, tf_index_down, stat_file_down, args.normfile)

    tock = time.time()-tick
    print('finished '+str(tock))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Identify transcription factors which may cause the different interactions between two Hi-C profiles.')

    parser.add_argument('-ci', '--controlindex', action='store', type=str, dest='c_index', required=True, help='HiC-pro index file for control sample', metavar='<str>')
    parser.add_argument('-cm', '--controlmatrix', action='store', type=str, dest='c_matrix', required=True, help='HiC-pro matrix file for control sample', metavar='<str>')
    parser.add_argument('-ti', '--treatindex', action='store', type=str, dest='t_index', required=True, help='HiC-pro index file for treatment sample', metavar='<str>')
    parser.add_argument('-tm', '--treatmatrix', action='store', type=str, dest='t_matrix', required=True, help='HiC-pro matrix file for treatment sample', metavar='<str>')

    parser.add_argument('-s', '--species',dest='species',type=str,metavar='<species>',choices=['hg38','mm10'],required=True,help='Species, please choose from "hg38" or "mm10".')

    parser.add_argument('-o', '--outdir', action='store', type=str, dest='outdir', help='output directory for Hi-C bart', metavar='<str>', default='BART3D_output/')
    parser.add_argument('-r', '--region', action='store', type=int, dest='region', help='Regions to expand when finding interactions', metavar='<int>', default=200000)

    args = parser.parse_args()
    if(len(sys.argv)) < 5: # two index, two matrix file, one species, are required
        parser.print_help()
        sys.exit(1)

    BART3D(args)
