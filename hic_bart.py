import os, sys
import argparse

import hicpro_matrix_to_juicer_matrix_format
import get_completed_viewregion_interaction
import append_all_bins_both_side_normalized_interaction
import hic_interaction_change

import utils

def main(args):
    os.makedirs(args.outdir, exist_ok=True)

    # get resolution of Hi-C data
    res_c = utils.get_resolution(args.c_index)
    res_t = utils.get_resolution(args.t_index)
    if res_c != res_t :
        sys.stderr.write("Error: Resolutions from two Hi-C index file are not the same!\n")
        sys.exit(1)
    args.resolution = res_c

    # TODO: according to different species, get chrom and length for the species
    # args.species


    # get prefix of file
    c_prefix = os.path.splitext(os.path.basename(args.c_index))[0]
    t_prefix = os.path.splitext(os.path.basename(args.t_index))[0]

    # output: {prefix}_{resolution}_{flag}_{chrom}.matrix
    sys.stdout.write("Step1: get original Hi-C matrix data from Juicer format..\n")
    hicpro_matrix_to_juicer_matrix_format.write_out_juicer_format_matrix(args.c_index, args.c_matrix, args.outdir, c_prefix, args.resolution, 'control')
    hicpro_matrix_to_juicer_matrix_format.write_out_juicer_format_matrix(args.t_index, args.t_matrix, args.outdir, t_prefix, args.resolution, 'treat')

    # output: {prefix}_{chrom}_res_{resolution}_view_region_{region}.csv
    sys.stdout.write("Step2: get complete Hi-C interaction within {} region.. \n".format(args.region))
    for chrom in utils.chroms:
        c_matrix_file = args.outdir+os.sep+'{}_{}_{}_{}.matrix'.format(c_prefix, args.resolution, 'control', chrom)
        c_matrix_df = utils.get_matrix_df(c_matrix_file, args.region)
        if c_matrix_df is None:
            continue
        get_completed_viewregion_interaction.write_out_interactions(c_matrix_df, args.region, args.outdir, c_prefix, args.resolution, chrom)
        
        t_matrix_file = args.outdir+os.sep+'{}_{}_{}_{}.matrix'.format(t_prefix, args.resolution, 'treat', chrom)
        t_matrix_df = utils.get_matrix_df(t_matrix_file, args.region)
        if t_matrix_df is None:
            continue
        get_completed_viewregion_interaction.write_out_interactions(t_matrix_df, args.region, args.outdir, t_prefix, args.resolution, chrom)

    
    # output: {prefix}_res{resolution}_view{region}_{chrom}.csv
    sys.stdout.write("Step3: mirror the interactions to +- regions/2.. \n")
    append_all_bins_both_side_normalized_interaction.write_out_binding_interactions_sep_chroms(c_prefix, args.region, args.resolution, args.outdir)
    append_all_bins_both_side_normalized_interaction.write_out_binding_interactions_sep_chroms(t_prefix, args.region, args.resolution, args.outdir)

    # output: {treat}_over_{control}_res{}_view{}.csv
    sys.stdout.write("Step4: pair test between control and treat..\n")
    hic_interaction_change.compare_hic_interaction([t_prefix, c_prefix], args.region, args.resolution, args.outdir, args.outdir)


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