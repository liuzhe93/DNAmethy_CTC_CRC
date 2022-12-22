setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/04_FunctionalEnrichmentAnalysis/02_CorrlationMethyExp")
rm(list = ls())

methy<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/Data_Collection/TCGA-COAD/methylation/TCGA_coad.met_uniq.csv",header = TRUE, row.names = 1)
dataexpr<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/Overlap_genes/FPKM_mRNAExp.Me_NoMe_selected.csv",header=TRUE,row.names = 1)
genes <- c("TRIP10", "NGFR", "SLC48A1", "SRMS")

methy_selected<-methy[genes,]
methy_selected_rmNA<-na.omit(methy_selected)
FPKM_selected<-dataexpr
FPKM_selected<-FPKM_selected[genes,]
dim(FPKM_selected)
#[1] 4 449
dim(methy_selected_rmNA)
#[1] 4 351
colnames(FPKM_selected)<-substr(colnames(FPKM_selected),1,16)
colnames(methy_selected_rmNA)<-substr(colnames(methy_selected_rmNA),1,16)
sameSamples<-intersect(colnames(FPKM_selected),colnames(methy_selected_rmNA))
length(sameSamples)
#[1] 261
FPKM_selected1=FPKM_selected[,sameSamples]
methy_selected1=methy_selected_rmNA[,sameSamples]
sameGenes<-intersect(rownames(FPKM_selected),rownames(methy_selected_rmNA))
length(sameGenes)
#[1] 4
FPKM_selected2=FPKM_selected1[sameGenes,]
methy_selected2=methy_selected1[sameGenes,]

library("ggExtra")
library("ggplot2")
library("ggpubr")
library("ggpmisc")
theme_set(ggpubr::theme_pubr()+
            theme(legend.position = "top"))
###############################################TRIP10############################################################
xval<-t(methy_selected2["TRIP10",])
yval<-t(FPKM_selected2["TRIP10",])
data_TRIP10<-cbind(xval,yval)
colnames(data_TRIP10)<-c("Methylation","Expression")
data_TRIP10<-as.data.frame(data_TRIP10)
pdf("TRIP10.pdf")
p<-ggscatter(data_TRIP10,x="Methylation",y="Expression",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
###############################################NGFR###############################################################
xval<-t(methy_selected2["NGFR",])
yval<-t(FPKM_selected2["NGFR",])
data_NGFR<-cbind(xval,yval)
colnames(data_NGFR)<-c("Methylation","Expression")
data_NGFR<-as.data.frame(data_NGFR)
pdf("NGFR.pdf")
p<-ggscatter(data_NGFR,x="Methylation",y="Expression",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
###############################################SLC48A1############################################################
xval<-t(methy_selected2["SLC48A1",])
yval<-t(FPKM_selected2["SLC48A1",])
data_SLC48A1<-cbind(xval,yval)
colnames(data_SLC48A1)<-c("Methylation","Expression")
data_SLC48A1<-as.data.frame(data_SLC48A1)
pdf("SLC48A1.pdf")
p<-ggscatter(data_SLC48A1,x="Methylation",y="Expression",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
###############################################SRMS###############################################################
xval<-t(methy_selected2["SRMS",])
yval<-t(FPKM_selected2["SRMS",])
data_SRMS<-cbind(xval,yval)
colnames(data_SRMS)<-c("Methylation","Expression")
data_SRMS<-as.data.frame(data_SRMS)
pdf("SRMS.pdf")
p<-ggscatter(data_SRMS,x="Methylation",y="Expression",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
