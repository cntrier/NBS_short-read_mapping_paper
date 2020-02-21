module load samtools
module load python3

### Calculate depth for each population bam file for 150bp library and all read lengths at MQ10 and MQ20


for pop in GWD FIN GIH CHS CLM; do
  samtools depth -Q 19 -b gene_only_problem_fixed.bed -a ../fasta_files/${pop}_150.bam > ${pop}_150_MQ20.txt
done


for pop in 70 100 150 250; do
   samtools depth -Q 19 -b gene_only_problem_fixed.bed -a ../fasta_files/${pop}.bam > ${pop}_MQ20.txt
done

for pop in GWD FIN GIH CHS CLM; do
  samtools depth -Q 9 -b gene_only_problem_fixed.bed -a ../fasta_files/${pop}_150.bam > ${pop}_150_MQ10.txt
done



for pop in 70 100 150 250; do
   samtools depth -Q 9 -b gene_only_problem_fixed.bed -a ../fasta_files/${pop}.bam > ${pop}_MQ10.txt
done
