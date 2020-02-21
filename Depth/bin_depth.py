### Bin depth calculations into regions
### Python3
### Usage python3 depth_calc.py <depth_input> <bed_file> <outfile_name>
### 12.3.2019

import pandas as pd
import sys
import math

## Import files

depth = pd.read_csv(sys.argv[1], sep='\t', header=0)
bed = pd.read_csv(sys.argv[2], sep='\t', header=None)
bed.columns=["chr", "start", "end", "gene"]


f= open(sys.argv[3],'w')

# Bin into regions

for index, row in bed.iterrows():
	chr = str(row[0])
	start = int(row[1])
	end = int(row[2])
	gene_name =str(row[3])
	bed_chr = depth[(depth['chr'] == chr)]
	print(bed_chr.head())
	gene = bed_chr[(bed_chr['pos'] > start) & (bed_chr['pos'] < end)]
	min=20
	print(gene)
	for index, bp in gene.iterrows():
		depth_ind = float(bp[3])
		chr_bp = str(bp[0])
		pos_bp = int(bp[1])
		gene_bp = gene_name 
		if depth_ind <= min :
			line="%s\t%i\t%d\t%s\n" % (chr_bp, pos_bp, depth_ind, gene_bp)
			f.writelines(line)
f.close()
