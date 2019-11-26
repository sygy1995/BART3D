import os,sys,argparse
import numpy as np
import pandas as pd
from scipy import stats

import BART3D.utils as utils

def compare_hic_interaction(control_np,treatment_np,resolution,file_up,file_down,chrom,species):
	compr_data_out_up = open(file_up,'a')
	compr_data_out_down = open(file_down,'a')

	bin_number = len(control_np)
	for i in range(bin_number):
	    stats_score,pvalue = stats.ttest_rel(treatment_np[i],control_np[i])
	    if np.isnan(pvalue):
	        pvalue = 1
	        stats_score = 0
	    start = i*resolution
	    compr_data_out_up.write('{}\t{}\t{}\t.\t{:.3f}\t.\n'.format(chrom, start, start+resolution, stats_score))
	    compr_data_out_down.write('{}\t{}\t{}\t.\t{:.3f}\t.\n'.format(chrom, start, start+resolution, -stats_score))
	sys.stdout.write("Written {} bins for chromosome {}..\n".format(i+1,chrom))