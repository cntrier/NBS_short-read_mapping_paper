#### Print average GEM score for exons of interest
#### Usage python umap_fix.py <gene>

import sys
import pandas as pd

#gene = sys.argv[1]
genes = ["ABCD1", "ACAD8", "ACADM", "ACADS", "ACADSB", "ACADVL", "ACAT1", "ADA", "ADK", "AHCY", "AK2", "ALDH7A1", "ARG1", "ARPC1B", "ARSA", "ASL", "ASS1", "ATM", "BCKDHA", "BCKDHB", "BCKDK", "BCL11B", "BLNK", "BTD", "BTK", "CBS", "CD247", "CD3D", "CD3E", "CD3G", "CD40LG", "CD79A", "CD79B", "CD8A", "CFTR", "CHD7", "CORO1A", "CPS1", "CPT1A", "CPT2", "CYP11A1", "CYP17A1", "CYP21A1P", "CYP21A2", "DBT", "DCLRE1C", "DKC1", "DLD", "DMD", "DNMT3B", "DOCK2", "DOCK8", "EFL1", "ETFB", "ETFDH", "EXTL3", "FAH", "FOXN1", "GAA", "GALE", "GALK1", "GALT", "GAMT", "GATA2", "GATM", "GCDH", "HADHA", "HADHB", "HBB", "HLCS", "HMGCL", "HMGCS2", "HPD", "IDUA", "IGHM", "IGLL1", "IKBKB", "IKZF1", "IL2RG", "IL7R", "IVD", "JAK3", "LAT", "LCK", "LIG4", "LMBRD1", "LRRC8A", "LRRC8B", "LRRC8C", "LRRC8D", "MAT1A", "MC2R", "MCCC1", "MCCC2", "MMAA", "MMAB", "MMACHC", "MMADHC", "MTHFR", "MTR", "MTRR", "MUT", "NADK2", "NAGS", "NBN", "NHEJ1", "NPC1", "NPC2", "NUDCD3", "OTC", "PAH", "PAX1", "PCBD1", "PCCA", "PCCB", "HPGD", "PGM3", "PIK3CD", "PIK3R1", "PNP", "POR", "PRKDC", "PSAP", "PSAT1", "PSPH", "PTPRC", "PTS", "QDPR", "RAC2", "RAG1", "RAG2", "RECQL4", "CARMIL2", "RMRP", "SBDS", "SLC16A1", "SLC22A5", "SLC25A13", "SLC25A15", "SLC25A20", "SLC46A1", "SLC52A2", "SLC52A3", "SLC6A8", "SMARCAL1", "SMN1", "SMN2", "STAT5B", "STIM1", "TAT", "TBX1", "TCF3", "TTC7A", "UNC119", "WAS", "ZAP70", "ZBTB24", "ETFA", "LRRC8E"]


for i in genes:
	gene = str(i)
	gem=pd.read_table("75kmer_GEM.bed", header=None)
	gem.columns = ["chr", "start", "stop","id", "score"]
	## Load exon data for each gene 
	exons = pd.read_table("%s_exons.bed" % gene, sep='\t')

	##### Bin GEM data into exons #### 

	## Break up to reduce file size 

	
	score = []

	gem['score'] = gem['score'].astype(float)

	for index, row in exons.iterrows():
		chr = str(row[0])
		start = float(row[1])
		end = float(row[2])
		gem_chr = gem[(gem['chr'].astype(str) == chr]
		y = gem_chr[(gem_chr['start'] < start ) & (gem_chr['stop'] > end)]
		if not y.empty:
			kscore = y['score'].iloc[0]
			score.append(kscore)
		else:
			x = gem_chr.groupby(pd.cut(gem_chr['start'], bins=[start, end])).score.mean()
			score.append(x[0])

	gem_score = pd.DataFrame(score)
	gem_score.columns= ['gem_score']
	gem_score=gem_score.dropna()

	if gem_score.empty:
		print(gene)
	else: 
		print(gem_score)
	### Put it all together in a final map quality file for each exon

	exons_map_gem = pd.concat([exons, gem_score], axis = 1)
	exons_map_gem.to_csv("%s_final_map_gem.txt" % gene, sep='\t', index= None)