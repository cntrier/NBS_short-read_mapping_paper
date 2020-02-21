### Collect info from camo regions and where they mapped to from different files
### Output a bed file for every RL
### Python3 
### 8/7/2019

import pandas as pd
import re
import pybedtools
import numpy as np

pops=['250','70', '100', '150']
for j in range(0,len(pops)):
	pop=pops[j]
	print(pop)
	### This file has information on where they incorrectly mapped to
	alts= pd.read_csv("%s_alternate_hits_50.txt" % pop, header=None, sep='\t')
	alts.columns= ("chr", "start", "end", "sim_chr", "sim_start", "sim_end", "map_chr", "map_start", "map_end")
	alts = alts.drop(["sim_chr", "sim_start", "sim_end"], axis = 1)
	alts = alts.drop_duplicates()

	### This file has the gene information for incorrect hits
	pseudo = pd.read_csv("gene_%s_pseudo_50.bed" % pop, header=None, sep='\t')
	pseudo.columns = ("hit_chr", "hit_start", "hit_end", "type", "gene")

	### This file has gene annotations for camoflagued regions
	gene_info = pd.read_csv("gene_%s_50.bed" % pop, header=0, sep='\t')
	gene_info= gene_info.drop(["gene_type"], axis=1)

	#### This file has exon annotations for camouflaged regions 
	exon_info = pd.read_csv("exons_%s_50.bed" % pop, header=None, sep='\t')
	exon_info.columns = ("chr", "start", "end", "exon_number")

	### Merge gene and exon info into one file
	annotations = pd.merge(gene_info, exon_info, how="left", on=(['chr','start', 'end']))
	
	print("####### Bed file annotation ####### ")

	annotations['exon_number'] = annotations['exon_number'].fillna(-1)
	annotations['exon_number'] = annotations['exon_number'].astype(int)
	annotations['exon_number']= annotations['exon_number'].astype(str)
	annotations['exon_number'] = annotations['exon_number'].replace('-1', 'Intronic')


	### Combine alternate hits with annotation information

	alts_annotated = pd.merge(alts, annotations, how="right", on=(['chr','start', 'end']))
	
	print("###### Merged alts and anotations #########")

	### Only deal with mappings to the primary assembly

	chr_names = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chr23", "chrX"]
	alts_annotated = alts_annotated[alts_annotated.map_chr.isin(chr_names)]

	### Combine pseudogene annotations with alternate hits file in duplicate and then drop the map_start stuff and then remove duplicate lines
	
	for index, row in pseudo.iterrows():
		hit_chr = str(row[0])
		hit_start = int(row[1])
		hit_end = int(row[2])
		hit_type = str(row[3])
		hit_name = str(row[4])
		bed_set=alts_annotated[(alts_annotated['map_chr']==hit_chr)]
		overlap = bed_set[(bed_set['map_start'] >= hit_start) & (bed_set['map_end'] <= hit_end)]
		overlap['hit_type'] = hit_type
		overlap['hit_name'] = hit_name
		overlap['hit_region'] = '{0}:{1}-{2}'.format(hit_chr, hit_start, hit_end)
		overlap=overlap.drop(["map_chr", "map_start", "map_end"], axis = 1)
		overlap= overlap.drop_duplicates()
		if index == 0:
			df= overlap
		elif index !=0:
			df= pd.concat([df, overlap], axis=0)

	print ('###### Hit Annotations Combined #####')
	print(df)

	#### Combine with information on total number of reads in that region

	count = pd.read_csv("%s_alternate_mapping_count.txt" % pop, header=None, sep='\t')
	count.columns=["chr", "start", "end", "count", "length", "rl"]
	df_all = df.merge(count, how="right", on=(['chr','start', 'end']))
	
	alt_count = pd.read_csv("%s_alternate_hit_count_10.txt" % pop, header=None, sep='\t')
	alt_count.columns=["chr", "start", "end", "alt_count", "length", "rl"]
	alt_count = alt_count.drop(["length", "rl"], axis = 1)

	final = df_all.merge(alt_count, how="right", on=(['chr','start', 'end']))
	final['percent_incorrect']= round(final['alt_count'] / final['count'] * 100, 2)
	print(final)
	print('#### Count File Added #####')
	
	regions_bed = final.drop(["hit_type", "hit_name", "hit_region", "length", "rl"], axis = 1)
	regions_bed = regions_bed.drop_duplicates()
	regions_bed['percent_incorrect']= regions_bed['alt_count'] / regions_bed['count'] * 100
	regions_bed=regions_bed.drop(["count", "alt_count"], axis =1)
	regions_bed.to_csv("%s_Annotated_Problem_Regions.txt" % pop, sep='\t', index = False)

	bytag = final.groupby(['chr','start', 'end', 'gene_name','exon_number','percent_incorrect'])
	print(bytag.head())
	
	f=open('%s_Alternate_Regions.bed' % pop, 'w')

	for name, group in bytag:
		
		print(name)
		
		chr = name[0]
		start = name[1]
		end = name[2]
		gene_name = str(name[3])
		exon = name[4]
		incorrect = name[5]

		hits = bytag.get_group(name)["hit_region"]
		hits = ','.join(hits)


		row = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(chr, start, end, gene_name, exon, incorrect, hits)
		print(row)
		f.writelines(row)
	f.close()	

	hit_annotations = final.drop(["length", "rl"], axis = 1)
	hit_annotations['percent_incorrect']= round(hit_annotations['alt_count'] / hit_annotations['count'] * 100,2)
	hit_annotations=hit_annotations.drop(["count", "alt_count"], axis =1)
	print(hit_annotations)
	hit_annotations.to_csv("%s_Annotated_Alternate_Regions.txt" % pop, sep='\t', index = False)