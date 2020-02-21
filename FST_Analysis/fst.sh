for pop in CHS FIN GWD CLM GIH; do
	./angsd/angsd -b ${pop}.list -anc GRCh38_full_analysis_set_plus_decoy_hla.fa -out ${pop} -dosaf 1 -gl 2
done	

./angsd/misc/realSFS CHS.saf.idx FIN.saf.idx -P 16 > CHS.FIN.ml
./angsd/misc/realSFS CHS.saf.idx GWD.saf.idx -P 16 > CHS.GWD.ml
./angsd/misc/realSFS CHS.saf.idx GIH.saf.idx -P 16 > CHS.GIH.ml
./angsd/misc/realSFS CHS.saf.idx CLM.saf.idx -P 16 > CHS.CLM.ml

./angsd/misc/realSFS GIH.saf.idx GWD.saf.idx -P 16 > GIH.GWD.ml
./angsd/misc/realSFS FIN.saf.idx GWD.saf.idx -P 16 > FIN.GWD.ml
./angsd/misc/realSFS GWD.saf.idx CLM.saf.idx -P 16 > GWD.CLM.ml

./angsd/misc/realSFS GIH.saf.idx CLM.saf.idx -P 16 > GIH.CLM.ml
./angsd/misc/realSFS GIH.saf.idx FIN.saf.idx -P 16 > GIH.FIN.ml

./angsd/misc/realSFS FIN.saf.idx CLM.saf.idx -P 16 > FIN.CLM.ml


./angsd/misc/realSFS fst index CHS.saf.idx FIN.saf.idx GWD.saf.idx -sfs CHS.FIN.ml -sfs CHS.GWD.ml -sfs FIN.GWD.ml -fstout CHS_FIN_GWD 
./angsd/misc/realSFS fst index CHS.saf.idx GIH.saf.idx CLM.saf.idx -sfs CHS.GIH.ml -sfs CHS.CLM.ml -sfs GIH.CLM.ml -fstout CHS_GIH_CLM
./angsd/misc/realSFS fst index GIH.saf.idx GWD.saf.idx CLM.saf.idx -sfs GIH.GWD.ml -sfs GIH.CLM.ml -sfs GWD.CLM.ml -fstout GIH_GWD_CLM
./angsd/misc/realSFS fst index GIH.saf.idx FIN.saf.idx CLM.saf.idx -sfs GIH.FIN.ml -sfs GIH.CLM.ml -sfs FIN.CLM.ml -fstout GIH_FIN_CLM

./angsd/misc/realSFS fst stats2 CHS_FIN_GWD.fst.idx -win 1000 -step 1000 -type 1 > CHS_FIN_GWD.fst.slide.txt
./angsd/misc/realSFS fst stats2 CHS_GIH_CLM.fst.idx -win 1000 -step 1000 -type 1 > CHS_GIH_CLM.fst.slide.txt
./angsd/misc/realSFS fst stats2 GIH_GWD_CLM.fst.idx -win 1000 -step 1000 -type 1 > GIH_GWD_CLM.fst.slide.txt
./angsd/misc/realSFS fst stats2 GIH_FIN_CLM.fst.idx -win 1000 -step 1000 -type 1 > GIH_FIN_CLM.fst.slide.txt

./angsd/misc/realSFS fst stats CHS_FIN_GWD.fst.idx > CHS_FIN_GWD_fst_global.txt
./angsd/misc/realSFS fst stats CHS_GIH_CLM.fst.idx > CHS_GIH_CLM_fst_global.txt
./angsd/misc/realSFS fst stats GIH_GWD_CLM.fst.idx > GIH_GWD_CLM_fst_global.txt
./angsd/misc/realSFS fst stats GIH_FIN_CLM.fst.idx > GIH_FIN_CLM_fst_global.txt