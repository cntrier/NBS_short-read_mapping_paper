module load samtools

### Map simulated reads genomes and combine into bam files for analysis

for ID in $(cat 5pop_sample_id.list); do
	for RL in 70 100 150 250; do
		bwa mem -M -t 16 -B 4 -O 6 -E 1 ../bam_files/GRCh38_full_analysis_set_plus_decoy_hla.fa ${ID}_${RL}_read1.fastq.gz ${ID}_${RL}_read2.fastq.gz | samtools sort | samtools view -1 - > ${ID}_${RL}.bam
	done
done



for ID in $(cat 5pop_sample_id.list); do
	samtools merge -f ${ID}.bam ${ID}_70.bam ${ID}_100.bam ${ID}_150.bam ${ID}_250.bam 
done


samtools merge -f CHS.bam \
HG00176.bam \
HG00419.bam \
HG00437.bam \
HG00531.bam \
HG00566.bam \
HG00608.bam \
HG00632.bam \
HG00654.bam \
HG00675.bam \
HG00428.bam 


samtools merge -f CLM.bam \
HG01119.bam \
HG01251.bam  \
HG01345.bam  \
HG01366.bam  \
HG01375.bam  \
HG01390.bam  \
HG01444.bam  \
HG01456.bam  \
HG01468.bam  \
HG01134.bam   


samtools merge -f FIN.bam \
HG00174.bam  \
HG00176.bam  \
HG00272.bam  \
HG00276.bam  \
HG00288.bam  \
HG00309.bam  \
HG00331.bam  \
HG00334.bam  \
HG00377.bam  \
HG00378.bam  


samtools merge -f GIH.bam \
NA20882.bam  \
NA20854.bam  \
NA20856.bam  \
NA20862.bam  \
NA20868.bam \
NA20874.bam \
NA20876.bam \
NA20881.bam \
NA20899.bam \
NA20902.bam


samtools merge -f GWD.bam \
HG02646.bam \
HG02667.bam \
HG02763.bam \
HG02772.bam \
HG02805.bam \
HG02837.bam \
HG02879.bam \
HG02891.bam \
HG03046.bam \
HG03539.bam

for RL in 70 100 150 250; do
	samtools merge -f ${RL}.bam \
	HG02646_${RL}.bam \
	HG02667_${RL}.bam \
	HG02763_${RL}.bam \
	HG02772_${RL}.bam \
	HG02805_${RL}.bam \
	HG02837_${RL}.bam \
	HG02879_${RL}.bam \
	HG02891_${RL}.bam \
	HG03046_${RL}.bam \
	HG03539_${RL}.bam \
	HG00176_${RL}.bam \
	HG00419_${RL}.bam \
	HG00437_${RL}.bam \
	HG00531_${RL}.bam \
	HG00566_${RL}.bam \
	HG00608_${RL}.bam \
	HG00632_${RL}.bam \
	HG00654_${RL}.bam \
	HG00675_${RL}.bam \
	HG00428_${RL}.bam \
	HG01119_${RL}.bam \
	HG01251_${RL}.bam \
	HG01345_${RL}.bam \
	HG01366_${RL}.bam \
	HG01375_${RL}.bam \
	HG01390_${RL}.bam \
	HG01444_${RL}.bam \
	HG01456_${RL}.bam \
	HG01468_${RL}.bam \
	HG00176_${RL}.bam \
	HG00272_${RL}.bam \
	HG00276_${RL}.bam \
	HG00288_${RL}.bam \
	HG00309_${RL}.bam \
	HG00331_${RL}.bam \
	HG00334_${RL}.bam \
	HG00377_${RL}.bam \
	HG00378_${RL}.bam \
	NA20882_${RL}.bam \
	NA20854_${RL}.bam \
	NA20856_${RL}.bam \
	NA20862_${RL}.bam \
	NA20868_${RL}.bam \
	NA20874_${RL}.bam \
	NA20876_${RL}.bam \
	NA20881_${RL}.bam \
	NA20899_${RL}.bam \
	NA20902_${RL}.bam 
done