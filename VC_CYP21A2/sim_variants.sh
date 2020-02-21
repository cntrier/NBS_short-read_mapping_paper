### Create pathogenic variants for CYP21A2 ###

module load bcftools/1.9 

### Read simulation with dwgsim

for i in {1..15}; do 
	echo $i
	bcftools consensus -I -s HG02763 -f HG02763.fa CYP21A2_${i}_hom.vcf.gz > CYP21A2_variant_${i}_hom.fa
	bcftools consensus -I -s HG02763 -f HG02763.fa CYP21A2_${i}_het.vcf.gz > CYP21A2_variant_${i}_het.fa
		
		
	dwgsim -1 150 -2 150 -d 350 -s 10 -C 40 -R 0 -r 0 -n 5 -F 0 -e 0.0024 -E 0.0024 -Q 5 -q a -X 0 -x ../fasta_files/2kb_problem_genes.bed CYP21A2_variant_${i}_hom.fa CYP21A2_variant${i}_150_hom_all
	dwgsim -1 150 -2 150 -d 350 -s 10 -C 40 -R 0 -r 0 -n 5 -F 0 -e 0.0024 -E 0.0024 -Q 5 -q a -X 0 -x ../fasta_files/2kb_problem_genes.bed CYP21A2_variant_${i}_het.fa CYP21A2_variant${i}_150_het_all
	dwgsim -1 150 -2 150 -d 555 -s 10 -C 40 -R 0 -r 0 -n 5 -F 0 -e 0.0024 -E 0.0024 -Q 5 -q a -X 0 -x ../fasta_files/2kb_problem_genes.bed CYP21A2_variant_${i}_hom.fa CYP21A2_variant${i}_150_hom_insert
    dwgsim -1 150 -2 150 -d 555 -s 10 -C 40 -R 0 -r 0 -n 5 -F 0 -e 0.0024 -E 0.0024 -Q 5 -q a -X 0 -x ../fasta_files/2kb_problem_genes.bed CYP21A2_variant_${i}_het.fa CYP21A2_variant${i}_150_het_insert

	gunzip CYP21A2_variant${i}_150_het_all.bwa.read1.fastq.gz
	gunzip CYP21A2_variant${i}_150_het_all.bwa.read2.fastq.gz

	gunzip CYP21A2_variant${i}_150_het_insert.bwa.read1.fastq.gz
 	gunzip CYP21A2_variant${i}_150_het_insert.bwa.read2.fastq.gz

    ALT=$(sed "${i}q;d" alt.txt)
	REF=$(sed "${i}q;d" ref.txt)
	
	echo $REF
	echo $ALT
	
	sed "s/N/${ALT}/g" CYP21A2_variant${i}_150_het_all.bwa.read1.fastq > CYP21A2_variant${i}_150_het_subbed.bwa.read1.fastq
	sed "s/N/${REF}/g" CYP21A2_variant${i}_150_het_all.bwa.read2.fastq > CYP21A2_variant${i}_150_het_subbed.bwa.read2.fastq
	
	sed "s/N/${ALT}/g" CYP21A2_variant${i}_150_het_insert.bwa.read1.fastq > CYP21A2_variant${i}_150_het_insert_subbed.bwa.read1.fastq
    sed "s/N/${REF}/g" CYP21A2_variant${i}_150_het_insert.bwa.read2.fastq > CYP21A2_variant${i}_150_het_insert_subbed.bwa.read2.fastq

	bgzip CYP21A2_variant${i}_150_het_subbed.bwa.read1.fastq
	bgzip CYP21A2_variant${i}_150_het_subbed.bwa.read2.fastq
	bgzip CYP21A2_variant${i}_150_het_insert_subbed.bwa.read1.fastq
    bgzip CYP21A2_variant${i}_150_het_insert_subbed.bwa.read2.fastq

done