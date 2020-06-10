
for i in {1..11}; do 

#bgzip ${i}.vcf
#tabix ${i}.vcf.gz

#bcftools consensus -I -s HG02763 -f HG02763_fixed.fa ${i}.vcf.gz > ${i}_hom.fa

../DWGSIM/dwgsim -1 150 -2 150 -d 350 -s 10 -C 40 -R 0 -r 0 -n 5 -F 0 -e 0.0024 -E 0.0024 -Q 5 -q a -X 0 -x 250_2kb_merged.bed ${i}_hom.fa ${i}_150_hom
../DWGSIM/dwgsim -1 150 -2 150 -d 555 -s 10 -C 40 -R 0 -r 0 -n 5 -F 0 -e 0.0024 -E 0.0024 -Q 5 -q a -X 0 -x 250_2kb_merged.bed ${i}_hom.fa ${i}_150_hom_insert

bwa mem -M -t 1 -B 4 -O 6 -E 1 -R "@RG\tID:VC\tSM:Seq1\tPL:ILLUMINA\tPI:15"  GRCh38_full_analysis_set_plus_decoy_hla.fa ${i}_150_hom.bwa.read1.fastq.gz ${i}_150_hom.bwa.read2.fastq.gz | samtools sort | samtools view -1 - > ${i}_150_hom.bam
bwa mem -M -t 1 -B 4 -O 6 -E 1 -R "@RG\tID:VC\tSM:Seq1\tPL:ILLUMINA\tPI:15"  GRCh38_full_analysis_set_plus_decoy_hla.fa ${i}_150_hom_insert.bwa.read1.fastq.gz ${i}_150_hom_insert.bwa.read2.fastq.gz  | samtools sort | samtools view -1 - > ${i}_150_hom_insert.bam

samtools index  ${i}_150_hom.bam
samtools index ${i}_150_hom_insert.bam
	
gatk HaplotypeCaller -R  GRCh38_full_analysis_set_plus_decoy_hla.fa -I ${i}_150_hom.bam -ERC GVCF -O ${i}_hom_.g.vcf
gatk HaplotypeCaller -R  GRCh38_full_analysis_set_plus_decoy_hla.fa -I ${i}_150_hom_insert.bam -ERC GVCF -O ${i}_hom_insert.g.vcf	

gatk GenotypeGVCFs -R  GRCh38_full_analysis_set_plus_decoy_hla.fa --variant ${i}_hom_.g.vcf -O ${i}_hom.vcf 
gatk GenotypeGVCFs -R  GRCh38_full_analysis_set_plus_decoy_hla.fa --variant ${i}_hom_insert.g.vcf -O ${i}_hom_insert.vcf

done

