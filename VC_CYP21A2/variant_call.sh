module load samtools
module load gatk/4.0
module load bwa/0.7.17

for i in {1..11}; do 
	bwa mem -M -t 1 -B 4 -O 6 -E 1 -R "@RG\tID:CYP21A2\tSM:Seq1\tPL:ILLUMINA\tPI:15" GRCh38_full_analysis_set_plus_decoy_hla.fa  CYP21A2_variant${i}_150_hom_all.bwa.read1.fastq.gz CYP21A2_variant${i}_150_hom_all.bwa.read2.fastq.gz | samtools sort | samtools view -1 - > CYP21A2_variant${i}_150_hom.bam
	bwa mem -M -t 1 -B 4 -O 6 -E 1 -R "@RG\tID:CYP21A2\tSM:Seq1\tPL:ILLUMINA\tPI:15" GRCh38_full_analysis_set_plus_decoy_hla.fa  CYP21A2_variant${i}_150_het_subbed.bwa.read1.fastq.gz CYP21A2_variant${i}_150_het_subbed.bwa.read2.fastq.gz  | samtools sort | samtools view -1 - > CYP21A2_variant${i}_150_het.bam
	
	samtools index CYP21A2_variant${i}_150_hom.bam
	samtools index CYP21A2_variant${i}_150_het.bam
	
	gatk HaplotypeCaller -R GRCh38_full_analysis_set_plus_decoy_hla.fa  -I CYP21A2_variant${i}_150_hom.bam -ERC GVCF -O CYP21A2_variant${i}_hom_noalt.g.vcf
	gatk HaplotypeCaller -R GRCh38_full_analysis_set_plus_decoy_hla.fa  -I CYP21A2_variant${i}_150_het.bam  -ERC GVCF -O CYP21A2_variant${i}_het.noalt.g.vcf
	
	gatk HaplotypeCaller --minimum-mapping-quality 10 -R GRCh38_full_analysis_set_plus_decoy_hla.fa  -I CYP21A2_variant${i}_150_hom.bam -ERC GVCF -O CYP21A2_variant${i}_hom_noalt_MQ10.g.vcf
	gatk HaplotypeCaller --minimum-mapping-quality 10 -R GRCh38_full_analysis_set_plus_decoy_hla.fa  -I CYP21A2_variant${i}_150_het.bam -ERC GVCF -O CYP21A2_variant${i}_het.noalt_MQ10.g.vcf

	gatk GenotypeGVCFs -R GRCh38_full_analysis_set_plus_decoy_hla.fa --variant CYP21A2_variant${i}_het.noalt.g.vcf -O CYP21A2_variant${i}_het.final_noalt_150.vcf 
	gatk GenotypeGVCFs -R GRCh38_full_analysis_set_plus_decoy_hla.fa --variant CYP21A2_variant${i}_hom_noalt.g.vcf -O CYP21A2_variant${i}_hom.final_noalt_150.vcf

	gatk GenotypeGVCFs -R GRCh38_full_analysis_set_plus_decoy_hla.fa --variant CYP21A2_variant${i}_het.noalt_MQ10.g.vcf -O CYP21A2_variant${i}_het.final_noalt_150_MQ10.vcf
	gatk GenotypeGVCFs -R GRCh38_full_analysis_set_plus_decoy_hla.fa --variant CYP21A2_variant${i}_hom_noalt_MQ10.g.vcf -O CYP21A2_variant${i}_hom.final_noalt_150_MQ10.vcf

	bwa mem -M -t 1 -B 4 -O 6 -E 1 -R "@RG\tID:CYP21A2\tSM:Seq1\tPL:ILLUMINA\tPI:15" GRCh38_full_analysis_set_plus_decoy_hla.fa  CYP21A2_variant${i}_150_hom_insert.bwa.read1.fastq.gz CYP21A2_variant${i}_150_hom_insert.bwa.read2.fastq.gz | samtools sort | samtools view -1 - > CYP21A2_variant${i}_150_hom_insert.bam
	bwa mem -M -t 1 -B 4 -O 6 -E 1 -R "@RG\tID:CYP21A2\tSM:Seq1\tPL:ILLUMINA\tPI:15" GRCh38_full_analysis_set_plus_decoy_hla.fa  CYP21A2_variant${i}_150_het_insert_subbed.bwa.read1.fastq.gz  CYP21A2_variant${i}_150_het_insert_subbed.bwa.read2.fastq.gz   | samtools sort | samtools view -1 - > CYP21A2_variant${i}_150_het_insert.bam

	samtools index CYP21A2_variant${i}_150_hom_insert.bam
	samtools index CYP21A2_variant${i}_150_het_insert.bam

	gatk HaplotypeCaller -R GRCh38_full_analysis_set_plus_decoy_hla.fa  -I CYP21A2_variant${i}_150_hom_insert.bam -ERC GVCF -O CYP21A2_variant${i}_150_hom_insert.g.vcf
	gatk HaplotypeCaller -R GRCh38_full_analysis_set_plus_decoy_hla.fa  -I CYP21A2_variant${i}_150_het_insert.bam -O CYP21A2_variant${i}_150_het_insert.g.vcf

	gatk GenotypeGVCFs -R GRCh38_full_analysis_set_plus_decoy_hla.fa --variant CYP21A2_variant${i}_150_het_insert.g.vcf -O CYP21A2_variant${i}_het.final_150_insert.vcf 
	gatk GenotypeGVCFs -R GRCh38_full_analysis_set_plus_decoy_hla.fa --variant CYP21A2_variant${i}_150_hom_insert.g.vcf -O CYP21A2_variant${i}_hom.final_150_insert.vcf
done