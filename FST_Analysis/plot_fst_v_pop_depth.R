### Plot Fst values across NBS genes

### Import libraries
library(ggplot2)
library("grid")
library("gridExtra")
library(plyr)
library(tidyverse)
library(ggpubr)
library("gridExtra")
library(RColorBrewer)
library(reshape2)

### Read in data

CHS_FIN_GWD = read.csv("CHS_FIN_GWD.fst.slide_150.txt", sep="\t", header=TRUE)
CHS_GIH_CLM = read.csv("CHS_GIH_CLM.fst.slide_150.txt", sep="\t", header=TRUE)
GIH_FIN_CLM = read.csv("GIH_FIN_CLM.fst.slide_150.txt", sep="\t", header=TRUE)
GIH_GWD_CLM = read.csv("GIH_GWD_CLM.fst.slide_150.txt", sep="\t", header=TRUE)

CHS_FIN_GWD = CHS_FIN_GWD %>% dplyr::select(chr, midPos, Fst01, Fst02, Fst12)  %>% mutate(region = paste(CHS_FIN_GWD$chr,":", CHS_FIN_GWD$midPos,sep="")) %>% 
  dplyr::rename(CHS_FIN=Fst01) %>% dplyr::rename(CHS_GWD=Fst02) %>% dplyr::rename(GWD_FIN=Fst12)
  
CHS_GIH_CLM = CHS_GIH_CLM  %>% dplyr::select(chr, midPos, Fst01, Fst02, Fst12)  %>% mutate(region = paste(CHS_GIH_CLM$chr,":", CHS_GIH_CLM$midPos,sep="")) %>% dplyr::rename(CHS_GIH=Fst01) %>% dplyr::rename(CHS_CLM=Fst02) %>% dplyr::rename(GIH_CLM=Fst12) 

GIH_FIN_CLM = GIH_FIN_CLM %>% dplyr::select(chr, midPos, Fst01, Fst12)  %>% mutate(region = paste(GIH_FIN_CLM$chr,":", GIH_FIN_CLM$midPos,sep="")) %>% dplyr::rename(GIH_FIN=Fst01) %>% dplyr::rename(FIN_CLM=Fst12)
  
GIH_GWD_CLM = GIH_GWD_CLM  %>% dplyr::select(chr, midPos, Fst01,Fst12)  %>% mutate(region = paste(GIH_GWD_CLM$chr,":", GIH_GWD_CLM$midPos,sep="")) %>% dplyr::rename(GIH_GWD=Fst01) %>% dplyr::rename(GWD_CLM=Fst12) 

fst= CHS_FIN_GWD %>% inner_join(CHS_GIH_CLM, by=c("region","chr", "midPos")) %>% distinct
fst  = fst %>% inner_join(GIH_FIN_CLM, by=c("region","chr", "midPos")) %>% distinct
fst  = fst %>% inner_join(GIH_GWD_CLM , by=c("region","chr", "midPos")) %>% distinct


write.table(fst, "Fst_pops_slide.txt", sep="\t", row.names = FALSE, quote=F)

## Depth files

FIN_depth= read.csv("../depth_2019_11_7_/FIN_depth_binned_150.txt", sep="\t", header = FALSE)
colnames(FIN_depth)= c("chr", "midPos", "FIN_depth")

CHS_depth= read.csv("../depth_2019_11_7_/CHS_depth_binned_150.txt", sep="\t", header = FALSE)
colnames(CHS_depth)= c("chr", "midPos", "CHS_depth")

GIH_depth= read.csv("../depth_2019_11_7_/GIH_depth_binned_150.txt", sep="\t", header = FALSE)
colnames(GIH_depth)= c("chr", "midPos", "GIH_depth")

CLM_depth= read.csv("../depth_2019_11_7_/CLM_depth_binned_150.txt", sep="\t", header = FALSE)
colnames(CLM_depth)= c("chr", "midPos", "CLM_depth")

GWD_depth= read.csv("../depth_2019_11_7_/GWD_depth_binned_150.txt", sep="\t", header = FALSE)
colnames(GWD_depth)= c("chr", "midPos", "GWD_depth")


GWD_depth = GWD_depth %>% dplyr::select(GWD_depth)
GIH_depth = GIH_depth %>% dplyr::select(GIH_depth)
CLM_depth = CLM_depth %>% dplyr::select(CLM_depth)
FIN_depth = FIN_depth %>% dplyr::select(FIN_depth)


depth = cbind(FIN_depth, CHS_depth, GIH_depth, CLM_depth, GWD_depth)

depth["diff_GWD_FIN"]=abs(depth$GWD_depth-depth$FIN_depth)
depth["diff_CHS_FIN"]=abs(depth$CHS_depth-depth$FIN_depth)
depth["diff_CHS_GWD"]=abs(depth$CHS_depth-depth$GWD_depth)
depth["diff_CHS_GIH"]=abs(depth$CHS_depth-depth$GIH_depth)
depth["diff_CHS_CLM"]=abs(depth$CHS_depth-depth$CLM_depth)
depth["diff_GIH_CLM"]=abs(depth$GIH_depth-depth$CLM_depth)
depth["diff_GIH_FIN"]=abs(depth$GIH_depth-depth$FIN_depth)
depth["diff_FIN_CLM"]=abs(depth$FIN_depth-depth$CLM_depth)
depth["diff_GIH_GWD"]=abs(depth$GIH_depth-depth$GWD_depth)
depth["diff_GWD_CLM"]=abs(depth$GWD_depth-depth$CLM_depth)
depth['midPos']= depth$midPos


all = fst %>% inner_join(depth, by=c("chr", "midPos")) %>% distinct


## Plot

col=brewer.pal(n = 5, name = "Set2")


CHS_FIN_plot = ggscatter(all, x = "CHS_FIN", y = "diff_CHS_FIN", 
                         color =  col[1], size = 0.5,
                         add = "reg.line", conf.int = TRUE,
                         add.params = list(color = "black",fill = "lightgray"), 
                         cor.coef = TRUE,
                         cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top"),
                         xlab = "FST", ylab = "Delta Coverage", title ="CHS_FIN") + coord_cartesian(ylim = c(0, 12), xlim= c(0, 0.8))

CHS_GWD_plot = ggscatter(all, x = "CHS_GWD", y = "diff_CHS_GWD", 
                         color =  col[1], size = 0.5,
                         add = "reg.line", conf.int = TRUE,
                         add.params = list(color = "black",fill = "lightgray"), 
                         cor.coef = TRUE,
                         cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top"),
                         xlab = "FST", ylab = "Delta Coverage", title ="CHS_GWD") + coord_cartesian(ylim = c(0, 12), xlim= c(0, 0.8))

GWD_FIN_plot = ggscatter(all, x = "GWD_FIN", y = "diff_GWD_FIN", 
                         color =  col[2], size = 0.5,
                         add = "reg.line", conf.int = TRUE,
                         add.params = list(color = "black",fill = "lightgray"), 
                         cor.coef = TRUE,
                         cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top"),
                         xlab = "FST", ylab = "Delta Coverage", title ="GWD_FIN") + coord_cartesian(ylim = c(0, 12), xlim= c(0, 0.8))

CHS_GIH_plot = ggscatter(all, x = "CHS_GIH", y = "diff_CHS_GIH", 
                         color =  col[2], size = 0.5,
                         add = "reg.line", conf.int = TRUE,
                         add.params = list(color = "black",fill = "lightgray"),  
                         cor.coef = TRUE,
                         cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top"),
                         xlab = "FST", ylab = "Delta Coverage", title ="CHS_GIH")+ coord_cartesian(ylim = c(0, 12), xlim= c(0, 0.8))

CHS_CLM_plot = ggscatter(all, x = "CHS_CLM", y = "diff_CHS_CLM", 
                         color =  col[3], size = 0.5,
                         add = "reg.line", conf.int = TRUE,
                         add.params = list(color = "black",fill = "lightgray"),  
                         cor.coef = TRUE,
                         cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top"),
                         xlab = "FST", ylab = "Delta Coverage", title ="CHS_CLM") + coord_cartesian(ylim = c(0, 12), xlim= c(0, 0.8))

GIH_CLM_plot = ggscatter(all, x = "GIH_CLM", y = "diff_GIH_CLM", 
                         color =  col[3], size = 0.5,
                         add = "reg.line", conf.int = TRUE,
                         add.params = list(color = "black",fill = "lightgray"), 
                         cor.coef = TRUE,
                         cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top"),
                         xlab = "FST", ylab = "Delta Coverage", title ="GIH_CLM") + coord_cartesian(ylim = c(0, 12), xlim= c(0, 0.8))

GIH_FIN_plot = ggscatter(all, x = "GIH_FIN", y = "diff_GIH_FIN", 
                         color =  col[4], size = 1,
                         add = "reg.line", conf.int = TRUE,
                         add.params = list(color = "black",fill = "lightgray"), 
                         cor.coef = TRUE,
                         cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top"),
                         xlab = "FST", ylab = "Delta Coverage", title ="GIH_FIN") + coord_cartesian(ylim = c(0, 12), xlim= c(0, 0.8))

FIN_CLM_plot = ggscatter(all, x = "FIN_CLM", y = "diff_FIN_CLM", 
                         color =  col[4], size = 0.5,
                         add = "reg.line", conf.int = TRUE,
                         add.params = list(color = "black",fill = "lightgray"), 
                         cor.coef = TRUE,
                         cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top"),
                         xlab = "FST", ylab = "Delta Coverage", title ="FIN_CLM") + coord_cartesian(ylim = c(0, 12), xlim= c(0, 0.8))

GIH_GWD_plot = ggscatter(all, x = "GIH_GWD", y = "diff_GIH_GWD", 
                         color =  col[5], size = 0.5,
                         add = "reg.line", conf.int = TRUE,
                         add.params = list(color = "black",fill = "lightgray"), 
                         cor.coef = TRUE,
                         cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top"),
                         xlab = "FST", ylab = "Delta Coverage", title ="GIH_GWD") + coord_cartesian(ylim = c(0, 12), xlim= c(0, 0.8))

GWD_CLM_plot = ggscatter(all, x = "GWD_CLM", y = "diff_GWD_CLM", 
                         color =  col[5], size = 0.5,
                         add = "reg.line", conf.int = TRUE,
                         add.params = list(color = "black",fill = "lightgray"),
                         cor.coef = TRUE,
                         cor.coeff.args = list(method = "spearman", label.x.npc = "left", label.y.npc = "top"),
                         xlab = "FST", ylab = "Delta Coverage", title ="GWD_CLM") + coord_cartesian(ylim = c(0, 12), xlim= c(0, 0.8))


grid.arrange(CHS_FIN_plot,CHS_GWD_plot,GWD_FIN_plot,CHS_GIH_plot,CHS_CLM_plot,GIH_CLM_plot, GIH_FIN_plot, FIN_CLM_plot, GIH_GWD_plot, GWD_CLM_plot, ncol = 2, nrow = 5)


### Put Results in table

CHS_avg=mean(all$CHS_depth, na.rm = TRUE)
CHS_sd=sd(all$CHS_depth, na.rm = TRUE)

GWD_avg=mean(all$GWD_depth, na.rm = TRUE)
GWD_sd=sd(all$GWD_depth, na.rm = TRUE)

GIH_avg=mean(all$GIH_depth, na.rm = TRUE)
GIH_sd=sd(all$GIH_depth, na.rm = TRUE)

CLM_avg=mean(all$CLM_depth, na.rm = TRUE)
CLM_sd=sd(all$CLM_depth, na.rm = TRUE)

FIN_avg=mean(all$FIN_depth, na.rm = TRUE)
FIN_sd=sd(all$FIN_depth, na.rm = TRUE)

depth_table=cbind(CHS_avg, CHS_sd, GWD_avg, GWD_sd, GIH_avg, GIH_sd, CLM_avg, CLM_sd, FIN_avg, FIN_sd)

write.table(depth_table, "Overall_Depth_by_Pop.txt", sep="\t", row.names = FALSE, quote = FALSE)



##Plot fst peaks plot
k <- all %>% filter(chr == "chr6" & midPos > 80106647	& midPos < 80346270) %>% ggplot(aes((midPos), GWD_FIN, colour = "red")) + geom_point(size=0.5)
k <- k + coord_cartesian(ylim = c(0, 1))  + facet_grid(~chr, scales = "free_x", space="free_x") + ylab("FST") + geom_abline(intercept = 0, slope = 0)

GWD_CLM_plot = ggscatter(all, x = "GWD_CLM", y = "diff_GWD_CLM", 
                         color =  col[5], size = 1,
                         add = "reg.line", conf.int = TRUE,
                         add.params = list(color = "black",fill = "lightgray"),  
                         cor.coef = TRUE, cor.method = "spearman",
                         xlab = "FST", ylab = "Delta Coverage", title ="GWD_CLM") + coord_cartesian(ylim = c(-10, 10))


### Plot exons

CHS_FIN_GWD = read.csv("CHS_FIN_GWD_exons.fst.slide.txt", sep="\t", header=TRUE)
CHS_GIH_CLM = read.csv("CHS_GIH_CLM_exons.fst.slide.txt", sep="\t", header=TRUE)
GIH_FIN_CLM = read.csv("GIH_FIN_CLM_exons.fst.slide.txt", sep="\t", header=TRUE)
GIH_GWD_CLM = read.csv("GIH_GWD_CLM_exons.fst.slide.txt", sep="\t", header=TRUE)

CHS_FIN_GWD = CHS_FIN_GWD %>% dplyr::select(chr, midPos, Fst01, Fst02, Fst12) %>% mutate(region = paste(CHS_FIN_GWD$chr,":", CHS_FIN_GWD$midPos,sep="")) %>% 
  dplyr::rename(CHS_FIN=Fst01) %>% dplyr::rename(CHS_GWD=Fst02) %>% dplyr::rename(GWD_FIN=Fst12) %>% distinct()

CHS_GIH_CLM = CHS_GIH_CLM %>% dplyr::select(chr, midPos, Fst01, Fst02, Fst12) %>% dplyr::rename(CHS_GIH=Fst01) %>% dplyr::rename(CHS_CLM=Fst02) %>% dplyr::rename(GIH_CLM=Fst12) %>% 
  mutate(region = paste(CHS_GIH_CLM$chr,":", CHS_GIH_CLM$midPos,sep="")) %>% dplyr::rename(midPos2=midPos) %>% dplyr::rename(chr2=chr)

GIH_FIN_CLM = GIH_FIN_CLM %>% dplyr::select(Fst01, Fst12)  %>% dplyr::rename(GIH_FIN=Fst01) %>% dplyr::rename(FIN_CLM=Fst12)  

GIH_GWD_CLM = GIH_GWD_CLM %>% dplyr::select(Fst01, Fst12) %>% dplyr::rename(GIH_GWD=Fst01) %>% dplyr::rename(GWD_CLM=Fst12)

df1 = cbind(CHS_GIH_CLM,GIH_FIN_CLM, GIH_GWD_CLM)

df = CHS_FIN_GWD %>% inner_join(df1, by="region") %>% distinct


k <- df %>% ggplot(aes((midPos), GWD_FIN, colour = "red")) + geom_point(size=0.5)
k <- k + coord_cartesian(ylim = c(0, 1))  + facet_grid(~chr, scales = "free_x", space="free_x") + ylab("FST") + geom_abline(intercept = 0, slope = 0)


