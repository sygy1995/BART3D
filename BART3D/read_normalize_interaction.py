import numpy as np
import pandas as pd
import bisect
import BART3D.utils as utils

def get_normalized_viewpoint_interaction(chrom,index,matrix_df,region,resolution,species):


	# save the start and end index of each chromosome into a dict
	if species=='hg38':
	    chroms = utils.chroms_hg38
	    chrom_sizes = utils.chrom_size_df_hg38
	elif species=='mm10':
	    chroms = utils.chroms_mm10
	    chrom_sizes = utils.chrom_size_df_mm10
	    
	chrome_index_range = [min(index[chrom].keys()),max(index[chrom].keys())]

	# for chrom in chroms:
	# assume that the matrix file is already sorted
	# bisect out the first and last id1 that is in this chromosome
	low1 = bisect.bisect_left(matrix_df['id1'],chrome_index_range[0])
	high1 = bisect.bisect_right(matrix_df['id1'],chrome_index_range[1])
	# save all interaction starts with this chromosome into a distinct dataframe
	matrix_df_intra = matrix_df[low1:high1]
	# sort this chromosomal specific interaction matrix by id2
	matrix_df_intra = matrix_df_intra.sort_values(['id2'])
	# reassign the index to be from 0 to length
	matrix_df_intra = matrix_df_intra.set_index(pd.Index(list(range(0,len(matrix_df_intra)))))
	# bisect out the first and last id2 that is in this chromosome
	low2 = bisect.bisect_left(matrix_df_intra['id2'],chrome_index_range[0])
	high2 = bisect.bisect_right(matrix_df_intra['id2'],chrome_index_range[1])
	# update this chromosomal specific interaction matrix with interaction that both id1 and id2 are in this chromosome
	matrix_df_intra = matrix_df_intra[low2:high2]
	# sort this chromosomal specific interaction matrix by id1 and then id2
	matrix_df_intra = matrix_df_intra.sort_values(by=['id1','id2'])
	# reassign the index to be from 0 to length
	matrix_df_intra = matrix_df_intra.set_index(pd.Index(list(range(0,len(matrix_df_intra)))))
	index_dict = index[chrom]


	# initialize the viewpoint interaction 
	df_cols = np.arange(-(region-resolution),region,resolution)
	chrom_size = chrom_sizes.at[chrom,'len']
	df_index = np.arange(0,chrom_size,resolution)
	#interaction_dataframe = pd.DataFrame(0,index=df_index,columns=df_cols)


	#interaction_list = list(repeat(list(repeat(0,len(df_cols))),len(df_index)))
	interaction_np = np.zeros((len(df_index),len(df_cols)))
	# write interactions into the dataframe

	matrix_list_intra = list(zip(matrix_df_intra['id1'],matrix_df_intra['id2'],matrix_df_intra['count']))
	zero_position = int(region/resolution)-1

	for line in matrix_list_intra:
	    id_a,id_b,score = int(line[0]),int(line[1]),float(line[2])
	    if index_dict[id_b]-index_dict[id_a] < region:
	        position_a = int(index_dict[id_a]/resolution)
	        position_b = int(index_dict[id_b]/resolution)
	#        interaction_list[position_a][position_b-position_a+zero_position] = score
	#        interaction_list[position_b][position_a-position_b+zero_position] = score
	        interaction_np[position_a,position_b-position_a+zero_position] = score
	        interaction_np[position_b,position_a-position_b+zero_position] = score
	interaction_normalized_np = interaction_np/(interaction_np.sum(0)/len(interaction_np))

	return interaction_normalized_np