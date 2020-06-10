#bedtools maskfasta -fi GRCh38_full_analysis_set_plus_decoy_hla.fa -bed alternate_hits_2kb_merged.bed -fo GRCh38_masked.fa

#bwa index -a bwtsw GRCh38_masked.fa

#samtools faidx GRCh38_masked.fa


#java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
#	REFERENCE=GRCh38_masked.fa \
#	OUTPUT=GRCh38_masked.fa.dict

for i in {1..11}; do 
	#../DWGSIM/dwgsim -1 150 -2 150 -d 350 -s 10 -C 40 -R 0 -r 0 -n 5 -F 0 -e 0.0024 -E 0.0024 -Q 5 -q a -X 0 -x low_and_alt_2kb_merged.bed ${i}_hom.fa ${i}_150_both	

	#bwa mem -M -t 1 -B 4 -O 6 -E 1 -R "@RG\tID:CYP21A2\tSM:Seq1\tPL:ILLUMINA\tPI:15" GRCh38_masked.fa ${i}_150_both.bwa.read1.fastq.gz ${i}_150_both.bwa.read2.fastq.gz | samtools sort | samtools view -1 - > ${i}_150_hom_masked.bam
	
	#samtools index ${i}_150_hom_masked.bam
	
	#gatk HaplotypeCaller -R GRCh38_masked.fa -I ${i}_150_hom_masked.bam -ERC GVCF --sample-ploidy 4 -O ${i}_150_both_masked_pl4.g.vcf

	gatk GenotypeGVCFs -R GRCh38_masked.fa --variant ${i}_150_both_masked_pl4.g.vcf -O ${i}_150_both_masked_pl4.vcf

done
