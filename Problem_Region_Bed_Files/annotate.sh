### Annotate alternate hits 

for POP in 70 100 150 250;  do
	cut -f 7,8,9 ${POP}_alternate_hits_50.txt | bedtools sort | bedtools merge -d 50  > ${POP}_pseudomap.bed
	cut -f 1,2,3 ${POP}_alternate_hits_50.txt | bedtools sort | bedtools merge -d 50 > ${POP}_camo_coords.bed
	sort-bed ${POP}_pseudomap.bed | bedmap --ec --echo --echo-map-id --delim '\t' ${POP}_pseudomap.bed hg38.p12_gene_sorted.bed | awk '$4!=""' | gsed 's/;/\t/g' | cut -f 1,2,3,6,7 | sed 's/gene_name=//g' > gene_${POP}_pseudo_50.bed
done

