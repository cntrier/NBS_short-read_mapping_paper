### Run ANGSD to calculate gentoype likelihoods for PCA analysis

angsd -GL 2 -out genolike -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam bam.list

angsd -b bam.list -nThreads 10 -out all -GL 2 -doMaf 2 -doMajorMinor 1 -doGeno 32 -doPost 1 -SNP_pval 1e-3 -nind 50 -P 8
