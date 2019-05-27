import sys,argparse
import os,glob
import numpy as np
import pandas as pd

import utils


def write_out_juicer_format_matrix(order_index_file,matrix_file,outdir,data,resolution,flag):

    # read index file
    with open(order_index_file) as order_index_inf: 
        index_df = pd.read_csv(order_index_inf,sep='\t',header=None)
    index_df.columns = ['chr','x','y','id']
    
    # split the original matrix by chroms
    for chrom in utils.chroms:
        outfile = open(outdir+os.sep+'{}_{}_{}_{}.matrix'.format(data,resolution,flag,chrom),'w')
        index_df_chr = index_df[index_df['chr']==chrom]
        index_df_chr.index = index_df_chr['id']
        index_dict = index_df_chr[['x']].to_dict()['x']
        keys = index_dict.keys()
        
        matrix_inf = open(matrix_file)
        line = matrix_inf.readline()
        while line:
            sline = line.strip().split()
            id_a,id_b,score = int(sline[0]),int(sline[1]),float(sline[2])
            if (id_a in keys) and (id_b in keys):
                outfile.write('{}\t{}\t{:.2f}\n'.format(index_dict[id_a],index_dict[id_b],score))
            line = matrix_inf.readline()
        matrix_inf.close()
        outfile.close()

def main(data,resolution):
    par_outdir = 'f1_juicer_formatted_hicpro_matrix'
    par_indir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/11_HiC_analysis/f1_preprocess/HiC_Pro/HiC_fastq_split'
    #data = 'A6010'
    #resolution = 5000
    outdir = '{}/{}/{}'.format(par_outdir,data,resolution)
    os.makedirs(outdir,exist_ok=True)
    index_file = par_indir+os.sep+'{}/{}_out/hic_results/matrix/fastq/raw/{}/fastq_{}_abs.bed'.format(data,data,resolution,resolution)
    raw_matrix = par_indir+os.sep+'{}/{}_out/hic_results/matrix/fastq/raw/{}/fastq_{}.matrix'.format(data,data,resolution,resolution)
    iced_matrix = par_indir+os.sep+'{}/{}_out/hic_results/matrix/fastq/iced/{}/fastq_{}_iced.matrix'.format(data,data,resolution,resolution)
    print(data,resolution,os.path.isfile(index_file),os.path.isfile(raw_matrix),os.path.isfile(iced_matrix))
    write_out_juicer_format_matrix(index_file,raw_matrix,outdir,data,resolution,'raw')
    #write_out_juicer_format_matrix(index_file,iced_matrix,outdir,data,resolution,'iced')

 

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--data', action = 'store', type = str,dest = 'data', help = 'input file of', metavar = '<int>')
    parser.add_argument('-r', '--resolution', action = 'store', type = int,dest = 'resolution', help = 'input file of', metavar = '<str>')
    #parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of', metavar = '<file>')
    #parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main(args.data,args.resolution)
