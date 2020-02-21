module load samtools
module load gatk/4.0
module load bwa/0.7.17
module load bedtools
module load picard-tools/2.17.6

bedtools maskfasta -fi GRCh38_full_analysis_set_plus_decoy_hla.fa -bed CYP_hits.bed -fo  GRCh38_primary_CYP21A1P-masked.fa
bwa index -a bwtsw  GRCh38_primary_CYP21A1P-masked.fa
samtools faidx  GRCh38_primary_CYP21A1P-masked.fa

java -jar /cluster/software/VERSIONS/picard-tools/2.17.6/libs/picard.jar CreateSequenceDictionary \
	REFERENCE= GRCh38_primary_CYP21A1P-masked.fa \ 
	OUTPUT= GRCh38_primary_CYP21A1P-masked.fa.dict

for i in {1..15}; do 
	/projects/cees/bin/DWGSIM/dwgsim -1 150 -2 150 -d 350 -s 10 -C 40 -R 0 -r 0 -n 5 -F 0 -e 0.0024 -E 0.0024 -Q 5 -q a -X 0 -x both.bed CYP21A2_variant_${i}_hom.fa CYP21A2_variant${i}_150_hom_both
	/projects/cees/bin/DWGSIM/dwgsim -1 150 -2 150 -d 350 -s 10 -C 40 -R 0 -r 0 -n 5 -F 0 -e 0.0024 -E 0.0024 -Q 5 -q a -X 0 -x both.bed CYP21A2_variant_${i}_het.fa CYP21A2_variant${i}_150_het_both

	gunzip CYP21A2_variant${i}_150_het_both.bwa.read1.fastq.gz
	gunzip CYP21A2_variant${i}_150_het_both.bwa.read2.fastq.gz

	bgzip CYP21A2_variant${i}_150_hom_both.bwa.read1.fastq
	bgzip CYP21A2_variant${i}_150_hom_both.bwa.read2.fastq
		
	ALT=$(sed "${i}q;d" alt.txt)
	REF=$(sed "${i}q;d" ref.txt)

	sed "s/N/${ALT}/g" CYP21A2_variant${i}_150_het_both.bwa.read1.fastq > CYP21A2_variant${i}_150_het_subbed_both.bwa.read1.fastq
	sed "s/N/${REF}/g" CYP21A2_variant${i}_150_het_both.bwa.read2.fastq > CYP21A2_variant${i}_150_het_subbed_both.bwa.read2.fastq

	bgzip CYP21A2_variant${i}_150_het_subbed_both.bwa.read1.fastq
	bgzip CYP21A2_variant${i}_150_het_subbed_both.bwa.read2.fastq

	bwa mem -M -t 1 -B 4 -O 6 -E 1 -R "@RG\tID:CYP21A2\tSM:Seq1\tPL:ILLUMINA\tPI:15"  GRCh38_primary_CYP21A1P-masked.fa CYP21A2_variant${i}_150_hom_both.bwa.read1.fastq.gz CYP21A2_variant${i}_150_hom_both.bwa.read2.fastq.gz | samtools sort | samtools view -1 - > CYP21A2_variant${i}_150_hom_masked.bam
	bwa mem -M -t 1 -B 4 -O 6 -E 1 -R "@RG\tID:CYP21A2\tSM:Seq1\tPL:ILLUMINA\tPI:15"  GRCh38_primary_CYP21A1P-masked.fa CYP21A2_variant${i}_150_het_subbed_both.bwa.read1.fastq.gz CYP21A2_variant${i}_150_het_subbed_both.bwa.read2.fastq.gz  | samtools sort | samtools view -1 - > CYP21A2_variant${i}_150_het_masked.bam
	
	samtools index CYP21A2_variant${i}_150_hom_masked.bam
	samtools index CYP21A2_variant${i}_150_het_masked.bam
	
	gatk HaplotypeCaller -R  GRCh38_primary_CYP21A1P-masked.fa -I CYP21A2_variant${i}_150_hom_masked.bam --sample-ploidy 4 -ERC GVCF -O CYP21A2_variant${i}_hom_noalt_pl2.g.vcf
	gatk HaplotypeCaller -R  GRCh38_primary_CYP21A1P-masked.fa -I CYP21A2_variant${i}_150_het_masked.bam --sample-ploidy 4 -ERC GVCF -O CYP21A2_variant${i}_het.noalt_pl2.g.vcf
	
	gatk HaplotypeCaller --minimum-mapping-quality 10 -R  GRCh38_primary_CYP21A1P-masked.fa -I CYP21A2_variant${i}_150_hom_masked.bam --sample_ploidy 4 -ERC GVCF -O CYP21A2_variant${i}_hom_noalt_pl2_MQ10.g.vcf
	gatk HaplotypeCaller --minimum-mapping-quality 10 -R  GRCh38_primary_CYP21A1P-masked.fa -I CYP21A2_variant${i}_150_het_masked.bam --sample_ploidy 4 -ERC GVCF -O CYP21A2_variant${i}_het.noalt_pl2_MQ10.g.vcf

	gatk GenotypeGVCFs -R  GRCh38_primary_CYP21A1P-masked.fa --variant CYP21A2_variant${i}_het.noalt_pl2.g.vcf -O CYP21A2_variant${i}_het.final_noalt_pl2_150.vcf 
	gatk GenotypeGVCFs -R  GRCh38_primary_CYP21A1P-masked.fa --variant CYP21A2_variant${i}_hom_noalt_pl2.g.vcf -O CYP21A2_variant${i}_hom.final_noalt_pl2_150.vcf

	gatk GenotypeGVCFs -R  GRCh38_primary_CYP21A1P-masked.fa --variant CYP21A2_variant${i}_het.noalt_pl2_MQ10.g.vcf -O CYP21A2_variant${i}_het.final_noalt_pl2_150_MQ10.vcf
	gatk GenotypeGVCFs -R  GRCh38_primary_CYP21A1P-masked.fa --variant CYP21A2_variant${i}_hom_noalt_pl2_MQ10.g.vcf -O CYP21A2_variant${i}_hom.final_noalt_pl2_150_MQ10.vcf

done