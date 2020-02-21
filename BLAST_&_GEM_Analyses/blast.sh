module load bedtools
module load blast+

## Take coordinate file of relevant genes for clinical diagnosis blast against the human genome

### Set up reference genome as repeat masked blast database

## Acquire dusty information

dustmasker -in GRCh38.p12.fa -infmt fasta \
-outfmt maskinfo_asn1_bin -out GRCh38.p12.asnb

## Create blast database with masking information

makeblastdb -in GRCh38.p12.fa -dbtype nucl \
  -mask_data GRCh38.p12.asnb -out GRCh38.p12_dusty

## Check dusty info has been applied to the database
blastdbcmd -db GRCh38.p12_dusty -info

for f in $(cat relevant_genes.txt); do
        bedtools getfasta -fi GRCh38.p12.fa -bed ${f}_exons.bed -fo ${f}_exons_buffer70.fa
        blastn -task blastn -db GRCh38.p12_dusty -query ${f}_exons_buffer70.fa -dust yes -num_threads 5 -out ${f}_blast_hits.txt -outfmt '6 qseqid sseqid qstart qend sstart send pident length mismatch evalue'
done

