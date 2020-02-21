module load R/3.5.2

gunzip all.geno.gz

ngsCovar -probfile all.geno -outfile pop.covar -nind 50 -nsites 100000 -call 0

Rscript -e 'write.table(cbind(seq(1,50),rep(1,50),c(rep("CHS",10),rep("CLM",10),rep("FIN",10),rep("GIH", 10),rep("GWD", 10))), row.names=F, sep="\t", col.names=c("FID","IID","CLUSTER"), file="pops.clst", quote=F)'

Rscript plotPCA.R -i pop.covar -c 1-2 -a pops.clst -o pop.pca.pdf
