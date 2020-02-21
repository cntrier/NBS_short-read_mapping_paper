### Get number of simulated reads from each problem regions
### Python3
### 8.2.2019


import pandas as pd
import re

pops=['70', '100', '150', '250']

MQ_threshold = 0

for j in range(0,len(pops)):
	pop=pops[j]

	locs= pd.read_csv("%s_sam_filtered_positions_10.txt" % pop, header=None, sep='\t')
	locs.columns=['sim_chr', 'sim_start', "sim_end"]

	bed=pd.read_csv("%s_camo_50.bed" % pop, sep="\t", header=None)
	bed.columns=['chr', 'start', 'end']

	f=open('%s_alternate_mapping_count.txt' % pop, 'w')
		
	for index, row in bed.iterrows():
		bed_chr=str(row[0])
		bed_start =int(row[1])
		bed_end =int(row[2])
		pos = range(bed_start, bed_end)
		length = bed_end - bed_start
		### Same chromosome
		bed_set=locs[(locs['sim_chr']==bed_chr)]
		bed_set_start = bed_set[bed_set['sim_start'].isin(pos)]
		bed_set_end = bed_set[bed_set['sim_end'].isin(pos)]
		inner = bed_set_start.merge(bed_set_end, on = ['sim_chr', 'sim_start', 'sim_end'], how = 'inner')
		outer = bed_set_start.merge(bed_set_end, on = ['sim_chr', 'sim_start', 'sim_end'], how = 'outer')
		print(inner)
		print(outer)
		final = pd.concat([inner, outer], axis = 0)

		if length < int(pop):
			bed_set_over = bed_set[((bed_set['sim_start'] < bed_start) & (bed_set['sim_end'] > bed_end))]
			final = pd.concat([final, bed_set_over], axis = 0)
	
		count = len(final.index)
		print(final)
		print(count)
		line = "%s\t%i\t%i\t%i\t%i\t%s\n" % (bed_chr, bed_start, bed_end, count, length, pop)
		print(line)
		f.writelines(line)
	f.close()