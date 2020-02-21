module load gatk
module load bcftools/1.9
module load samtools

### Use VCF file from 1000 genomes to create fasta files per individual 

for j in $(cat pops.txt); do 
	java -cp GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
    -R GRCh38_full_analysis_set_plus_decoy_hla.fa \
    -V ${j}_chr1.vcf \
    -V ${j}_chr2.vcf \
    -V ${j}_chr3.vcf \
    -V ${j}_chr4.vcf \
    -V ${j}_chr5.vcf \
    -V ${j}_chr6.vcf \
    -V ${j}_chr7.vcf \
    -V ${j}_chr8.vcf \
    -V ${j}_chr9.vcf \
    -V ${j}_chr10.vcf \
	-V ${j}_chr11.vcf \
	-V ${j}_chr12.vcf \
	-V ${j}_chr13.vcf \
	-V ${j}_chr14.vcf \
	-V ${j}_chr15.vcf \
	-V ${j}_chr16.vcf \
	-V ${j}_chr17.vcf \
	-V ${j}_chr18.vcf \
	-V ${j}_chr19.vcf \
	-V ${j}_chr20.vcf \
	-V ${j}_chr21.vcf \
	-V ${j}_chr22.vcf \
	-V ${j}_chrX.vcf \
    -out ${j}_all.vcf \
    -assumeSorted

	bgzip ${j}_all.vcf
	tabix ${j}_all.vcf.gz
	
	for k in $(cat ${j}.pop.txt); do
		bcftools consensus -I -s ${k} -f GRCh38_full_analysis_set_plus_decoy_hla.fa ${j}_all.vcf.gz > ${k}.fa
	done
done