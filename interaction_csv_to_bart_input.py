import sys,argparse
import os,glob
import numpy as np
import pandas as pd
from scipy import stats
import re,bisect

import utils

def main(viewregion,normalization,chrom):
    indir = 'f1_UDHS_interaction_change'
    outdir = 'f2_UDHS_bart_input'
    os.makedirs(outdir,exist_ok=True) 
    
    infiles = glob.glob(indir+'/*.csv')#;print(infiles);exit()
    for infile in infiles:
        basename = os.path.basename(infile).split('.csv')[0]
        with open(infile) as inf:
            df = pd.read_csv(inf,sep='\t',index_col=0)
        df = df[['stats']]
        df = df.fillna(0)
        df.to_csv(outdir+os.sep+'{}_Stats_Interaction_IncreasedTop.txt'.format(basename),header=None,sep='\t')
        # use the opposite values of stats scores
        df['stats'] = -1*df['stats']
        df.to_csv(outdir+os.sep+'{}_NegStats_Interaction_DecreasedTop.txt'.format(basename),header=None,sep='\t')
        #print(df);exit()



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
