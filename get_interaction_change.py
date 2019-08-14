import os,sys,argparse
import numpy as np
import pandas as pd
from scipy import stats

import utils

def compare_hic_interaction(control_np,treatment_np,resolution,file,chrom,species):
	compr_data_out = open(file,'a')

	bin_number = len(control_np)
	for i in range(bin_number):
	    stats_score,pvalue = stats.ttest_rel(treatment_np[i],control_np[i])
	    if np.isnan(pvalue):
	        pvalue = 1
	    start = i*resolution
	    compr_data_out.write('{}\t{}\t{}\t.\t{:.3f}\t.\n'.format(chrom, start, start+resolution, -np.log10(pvalue)))
	sys.stdout.write("Written {} bins for chromosome {}..\n".format(i+1,chrom))