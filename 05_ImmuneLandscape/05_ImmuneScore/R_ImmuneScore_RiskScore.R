setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/05_ImmuneLandscape/05_ImmuneScore")
rm(list = ls())

library(utils)
rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
#help(package="estimate")

risk<-read.table("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/04_multivariate_Cox_regression_analysis/risk.txt",
                 sep="\t", header=T, row.names = 1)
head(risk)
risk$sample_name<-row.names(risk)
risk$sample_name<-substr(risk$sample_name,1,12)
risk<-subset(risk, select = c(riskScore, risk, sample_name))

dataexpr<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/Overlap_genes/FPKM_mRNAExp.Me_NoMe_selected.csv",
                   header=TRUE, row.names = 1)
write.table(dataexpr, "dataexp.txt", sep = "\t", quote = F)

filterCommonGenes(input.f="dataexp.txt", 
                  output.f="COAD_59427genes.gct", 
                  id="GeneSymbol")
estimateScore(input.ds = "COAD_59427genes.gct",
              output.ds="COAD_estimate_score.gct", 
              platform="illumina")
scores=read.table("COAD_estimate_score.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
dim(scores)
#[1] 449   3
scores<-as.data.frame(scores)
scores$sample_name<-row.names(scores)
scores$sample_name<-substr(scores$sample_name,1,12)
for(i in 1:449){
  scores$sample_name[i]<-gsub("\\.", "-", scores$sample_name[i])
}
head(scores)
#                            StromalScore ImmuneScore ESTIMATEScore  sample_name
#TCGA.3L.AA1B.01A.11R.A37K.07   -448.06076    275.9807      -172.080 TCGA-3L-AA1B
#TCGA.4N.A93T.01A.11R.A37K.07  -1626.44660   -305.2113     -1931.658 TCGA-4N-A93T
#TCGA.5M.AATE.01A.11R.A41B.07  -1001.86366   -138.5890     -1140.453 TCGA-5M-AATE
#TCGA.A6.2672.01B.03R.2302.07     10.72671   1278.8947      1289.621 TCGA-A6-2672
#TCGA.A6.2672.01A.01R.0826.07    177.24664   1854.3007      2031.547 TCGA-A6-2672
#TCGA.A6.2676.01A.01R.0826.07    -50.55470   1578.7779      1528.223 TCGA-A6-2676

merged_results<-merge(risk, scores, by = "sample_name")
dim(merged_results)
#[1] 577   6
write.csv(merged_results, "estimateScore.csv", quote = F)
#merged_results<-read.csv("estimateScore.csv", header = T, row.names = 1)

data1<-subset(merged_results, select = c("risk", "riskScore", "StromalScore", "ImmuneScore", "ESTIMATEScore"))
data1_high<-subset(data1,risk == "high")
data1_low<-subset(data1,risk == "low")
library("ggExtra")
library("ggplot2")
library("ggpubr")
library("ggpmisc")
theme_set(ggpubr::theme_pubr()+
            theme(legend.position = "top"))
xval<-data1_high$riskScore
yval<-data1_high$ImmuneScore
data_riskScore_ImmuneScore<-cbind(xval,yval)
colnames(data_riskScore_ImmuneScore)<-c("RiskScore","ImmuneScore")
data_riskScore_ImmuneScore<-as.data.frame(data_riskScore_ImmuneScore)
pdf("Cor_RiskScore_ImmuneScore_HighGroup.pdf")
p<-ggscatter(data_riskScore_ImmuneScore,x="RiskScore",y="ImmuneScore",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()

xval<-data1_low$riskScore
yval<-data1_low$ImmuneScore
data_riskScore_ImmuneScore<-cbind(xval,yval)
colnames(data_riskScore_ImmuneScore)<-c("RiskScore","ImmuneScore")
data_riskScore_ImmuneScore<-as.data.frame(data_riskScore_ImmuneScore)
pdf("Cor_RiskScore_ImmuneScore_LowGroup.pdf")
p<-ggscatter(data_riskScore_ImmuneScore,x="RiskScore",y="ImmuneScore",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()


library(tinyarray)
library(tidyverse)
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
library(ggpubr)

for(i in 1:nrow(merged_results)){
  merged_results$sample_name[i] <- paste(merged_results$sample_name[i], i, sep = "_")
}
rownames(merged_results)<-merged_results$sample_name
merged_results<-merged_results[,-1]

dat <- merged_results %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Score,-Sample)

dat$sample_name<-dat$Sample
dat$sample_name<-substr(dat$sample_name,1,12)
merged_dat<-merge(dat, risk, by = "sample_name")
merged_dat$Group = merged_dat$risk
merged_dat$Group<-ifelse(merged_dat$Group=="high","2_high","1_low")
############################StromalScore#################################################
merged_dat_StromalScore <- subset(merged_dat, Cell_type == "StromalScore")
merged_dat_StromalScore <- subset(merged_dat_StromalScore, select = c("sample_name","Score","Group"))
merged_dat_StromalScore$Score <- as.numeric(merged_dat_StromalScore$Score)
pdf("Boxplot_StromalScore.pdf")
ggplot(merged_dat_StromalScore,aes(Group,Score,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Risk Group", y = "Stromal Score") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test")
dev.off()
############################ImmuneScore#################################################
merged_dat_ImmuneScore <- subset(merged_dat, Cell_type == "ImmuneScore")
merged_dat_ImmuneScore <- subset(merged_dat_ImmuneScore, select = c("sample_name","Score","Group"))
merged_dat_ImmuneScore$Score <- as.numeric(merged_dat_ImmuneScore$Score)
pdf("Boxplot_ImmuneScore.pdf")
ggplot(merged_dat_ImmuneScore,aes(Group,Score,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Risk Group", y = "Immune Score") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test")
dev.off()
############################ESTIMATEScore#################################################
merged_dat_ESTIMATEScore <- subset(merged_dat, Cell_type == "ESTIMATEScore")
merged_dat_ESTIMATEScore <- subset(merged_dat_ESTIMATEScore, select = c("sample_name","Score","Group"))
merged_dat_ESTIMATEScore$Score <- as.numeric(merged_dat_ESTIMATEScore$Score)
pdf("Boxplot_ESTIMATEScore.pdf")
ggplot(merged_dat_ESTIMATEScore,aes(Group,Score,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Risk Group", y = "ESTIMATE Score") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test")
dev.off()
################################riskScore#################################################
merged_dat_riskScore <- subset(merged_dat, Cell_type == "riskScore")
merged_dat_riskScore <- subset(merged_dat_riskScore, select = c("sample_name","Score","Group"))
merged_dat_riskScore$Score <- as.numeric(merged_dat_riskScore$Score)
pdf("Boxplot_riskScore.pdf")
ggplot(merged_dat_riskScore,aes(Group,Score,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Risk Group", y = "risk Score") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test")
dev.off()

