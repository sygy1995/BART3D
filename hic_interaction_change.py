import os,sys,argparse
import numpy as np
import pandas as pd
from scipy import stats

import utils

# def compare_hic_interaction(compr_data,viewregion,resolution,hic_normalized_interaction_info_dir,outdir):
#     treat,ctrl = compr_data[0],compr_data[1]
#     compr_data_out = open(outdir+os.sep+'{}_over_{}_res{}_view{}.bed'.format(treat,ctrl,resolution,viewregion),'w')

#     # output data: chr start end . -10log(pvalue) .
#     # compr_data_out.write('{}\t{}\t{}\t{}\t{}\n'.format('id','stats','pvalue','treat_mean','ctrl_mean'))

#     columns = np.arange(-1*viewregion+resolution,viewregion,resolution)       
    
#     for chrom in utils.chroms:
#         treat_hic_file = hic_normalized_interaction_info_dir+os.sep+'{}_res{}_view{}_{}.csv'.format(treat,resolution,viewregion,chrom)
#         ctrl_hic_file = hic_normalized_interaction_info_dir+os.sep+'{}_res{}_view{}_{}.csv'.format(ctrl,resolution,viewregion,chrom)
        
#         if not (os.path.isfile(treat_hic_file) and os.path.isfile(ctrl_hic_file)):
#             continue

#         with open(treat_hic_file) as treat_inf, open(ctrl_hic_file) as ctrl_inf:
#             treat_df = pd.read_csv(treat_inf,index_col=0)
#             ctrl_df = pd.read_csv(ctrl_inf,index_col=0)
        
#         #####
#         # get paired t-test pvalue
#         #####
#         if treat_df.index.tolist() != ctrl_df.index.tolist():
#             sys.stderr.write("Error: hic interaction has different index between normal and treat! \n")
#             sys.exit(1)

#         for index in treat_df.index.tolist():
#             view_pos_treat_data = np.nan_to_num(treat_df.loc[index][map(str,columns)].values)
#             view_pos_ctrl_data = np.nan_to_num(ctrl_df.loc[index][map(str,columns)].values)
#             stats_score,pvalue = stats.ttest_rel(view_pos_treat_data,view_pos_ctrl_data)
#             if np.isnan(pvalue):
#                 pvalue = 1
#             treat_mean = np.mean(view_pos_treat_data)
#             ctrl_mean = np.mean(view_pos_ctrl_data)
#             # output data: chr start end . -10log(pvalue) .
#             compr_data_out.write('{}\t{}\t{}\t.\t{:.3f}\t.\n'.format(chrom, index, index+resolution, -np.log10(pvalue)))
#             # compr_data_out.write('{}\t{:.3f}\t{:.3e}\t{:.3f}\t{:.3f}\n'.format(udhs_id,stats_score,pvalue,treat_mean,ctrl_mean))
#     compr_data_out.close()

def compare_hic_interaction(compr_data,viewregion,resolution,hic_normalized_interaction_info_dir,outdir):
    treat,ctrl = compr_data[0],compr_data[1]
    compr_data_out = open(outdir+os.sep+'{}_over_{}_res{}_view{}.bed'.format(treat,ctrl,resolution,viewregion),'w')

    # output data: chr start end . -10log(pvalue) .
    # compr_data_out.write('{}\t{}\t{}\t{}\t{}\n'.format('id','stats','pvalue','treat_mean','ctrl_mean'))

    columns = np.arange(-1*viewregion+resolution,viewregion,resolution)       
    print(utils.chroms)
    for chrom in utils.chroms:
        treat_hic_file = hic_normalized_interaction_info_dir+os.sep+'{}_res{}_view{}_{}.csv'.format(treat,resolution,viewregion,chrom)
        ctrl_hic_file = hic_normalized_interaction_info_dir+os.sep+'{}_res{}_view{}_{}.csv'.format(ctrl,resolution,viewregion,chrom)
        
        if not (os.path.isfile(treat_hic_file) and os.path.isfile(ctrl_hic_file)):
            continue

        with open(treat_hic_file) as treat_inf, open(ctrl_hic_file) as ctrl_inf:
            treat_lines = treat_inf.readlines()
            ctrl_lines = ctrl_inf.readlines()
            n = len(treat_lines)
            for i in range(1,n):
                treat_list = list(map(float,treat_lines[i].strip().split(',')))
                ctrl_list = list(map(float,ctrl_lines[i].strip().split(',')))
                index = treat_list[0]
                view_pos_treat_data = treat_list[1:]
                view_pos_ctrl_data = ctrl_list[1:]
                stats_score,pvalue = stats.ttest_rel(view_pos_treat_data,view_pos_ctrl_data)
                if np.isnan(pvalue):
                    pvalue = 1
                    # output data: chr start end . -10log(pvalue) .
                compr_data_out.write('{}\t{}\t{}\t.\t{:.3f}\t.\n'.format(chrom, index, index+resolution, -np.log10(pvalue)))
    compr_data_out.close()
            
def main(viewregion,normalization,chrom):
    outdir = 'f1_UDHS_interaction_change'
    os.makedirs(outdir,exist_ok=True) 
    
    hic_normalized_interaction_info_dir = '/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/11_HiC_analysis/f1_preprocess/union_binding_processed_data/panos/union_binding_view2M_interaction/f1_append_2M_BS_normalized_interaction_sep_chr'
    udhs_file='/nv/vol190/zanglab/zw5j/data/unionDHS/hg38_unionDHS_fc4_50merge.bed'
    with open(udhs_file) as udhs_inf:
        udhs_df = pd.read_csv(udhs_inf,header=None,index_col=3,sep='\t')
    udhs_df.columns = ['chr','start','end','len','strand']
    
    for compr_data in [['Jurkat','A6010'],['Cutll1','A6010']]:
        for resolution in [5000]:          
            compare_hic_interaction(udhs_df,compr_data,normalization,viewregion,resolution,hic_normalized_interaction_info_dir,outdir)



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
