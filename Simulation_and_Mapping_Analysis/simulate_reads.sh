module load samtools
module load bwa/0.7.17


for ID in $(cat 5pop_sample_id.list); do
	sed '/^>/ s/ .*//' ${k}.fa > ${k}_fixed.fa
	sed 's/R/A/g' ${ID}_fixed.fa | sed 's/Y/C/g' | sed 's/S/G/g' | sed 's/W/A/g' | sed 's/K/G/g' | sed 's/M/A/g' > ${ID}_fixed_1.fa
	sed 's/R/G/g' ${ID}_fixed.fa | sed 's/Y/T/g' | sed 's/S/C/g' | sed 's/W/T/g' | sed 's/K/T/g' | sed 's/M/C/g' > ${ID}_fixed_2.fa

	dwgsim -1 100 -2 100 -d 250 -s 10 -C 20 -R 0 -r 0 -n 5 -F 0 -e 0.0024 -E 0.0024 -Q 5 -q a -X 0 -x 2kb_problem_genes.bed ${ID}_fixed_1.fa ${ID}_100_1
	dwgsim -1 70 -2 70 -d 190 -s 10 -C 20 -R 0 -r 0 -n 5 -F 0 -e 0.0024 -E 0.0024 -Q 5 -q a -X 0 -x 2kb_problem_genes.bed ${ID}_fixed_1.fa ${ID}_70_1
	dwgsim -1 150 -2 150 -d 350 -s 10 -C 20 -R 0 -r 0 -n 5 -F 0 -e 0.0024 -E 0.0024 -Q 5 -q a -X 0 -x 2kb_problem_genes.bed ${ID}_fixed_1.fa ${ID}_150_1
	dwgsim -1 250 -2 250 -d 550 -s 10 -C 20 -R 0 -r 0 -n 5 -F 0 -e 0.0024 -E 0.0024 -Q 5 -q a -X 0 -x 2kb_problem_genes.bed ${ID}_fixed_1.fa ${ID}_250_1
	
	dwgsim -1 100 -2 100 -d 250 -s 10 -C 20 -R 0 -r 0 -n 5 -F 0 -e 0.0024 -E 0.0024 -Q 5 -q a -X 0 -x 2kb_problem_genes.bed ${ID}_fixed_2.fa ${ID}_100_2
	dwgsim -1 70 -2 70 -d 190 -s 10 -C 20 -R 0 -r 0 -n 5 -F 0 -e 0.0024 -E 0.0024 -Q 5 -q a -X 0 -x 2kb_problem_genes.bed ${ID}_fixed_2.fa ${ID}_70_2
	dwgsim -1 150 -2 150 -d 350 -s 10 -C 20 -R 0 -r 0 -n 5 -F 0 -e 0.0024 -E 0.0024 -Q 5 -q a -X 0 -x 2kb_problem_genes.bed ${ID}_fixed_2.fa ${ID}_150_2
	dwgsim -1 250 -2 250 -d 550 -s 10 -C 20 -R 0 -r 0 -n 5 -F 0 -e 0.0024 -E 0.0024 -Q 5 -q a -X 0 -x 2kb_problem_genes.bed ${ID}_fixed_2.fa ${ID}_250_2

done

for ID in $(cat 5pop_sample_id.list); do
	for RL in 70 100 150 250; do
		
		gunzip ${ID}_${RL}_1.bwa.read1.fastq.gz
		gunzip ${ID}_${RL}_1.bwa.read2.fastq.gz

		gunzip ${ID}_${RL}_2.bwa.read1.fastq.gz
		gunzip ${ID}_${RL}_2.bwa.read2.fastq.gz

		cat ${ID}_${RL}_1.bwa.read1.fastq ${ID}_${RL}_2.bwa.read1.fastq > ${ID}_${RL}_read1.fastq
		cat ${ID}_${RL}_1.bwa.read2.fastq ${ID}_${RL}_2.bwa.read2.fastq > ${ID}_${RL}_read2.fastq

		bgzip ${ID}_${RL}_read1.fastq
		bgzip ${ID}_${RL}_read2.fastq

	done
done