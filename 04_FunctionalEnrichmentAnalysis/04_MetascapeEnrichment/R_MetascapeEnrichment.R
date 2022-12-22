setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/04_FunctionalEnrichmentAnalysis/04_MetascapeEnrichment")
rm(list = ls())

risk<-read.table("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/04_multivariate_Cox_regression_analysis/risk.txt",
                 sep="\t",header=T, row.names = 1)
head(risk)
risk$sample_name<-row.names(risk)
risk$sample_name<-substr(risk$sample_name,1,12)
risk<-subset(risk, select = c(risk, sample_name))

dataexpr<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/Metastatic_nonMetastatic/dataexpr.csv",header=T,row.names=1)
dataexpr_t<-t(dataexpr)
dataexpr_t<-as.data.frame(dataexpr_t)
dataexpr_t$sample_name<-row.names(dataexpr_t)
dataexpr_t$sample_name<-substr(dataexpr_t$sample_name,1,12)
dataexpr_t$sample_name<-gsub("\\.","-",dataexpr_t$sample_name)
merged_data<-merge(dataexpr_t, risk, by = "sample_name")
dim(dataexpr_t)
#[1]   449 59428
dim(risk)
#[1] 521   2
dim(merged_data)
#[1]   577 59429

merged_data<-merged_data[order(merged_data$risk),]
merged_data[1:6,59426:59429]
#  ZYXP1 ZZEF1 ZZZ3 risk
#1     0  3107 2105 high
#2     0  2782  928 high
#3     0  2741 1582 high
#4     0  2274 1311 high
#5     0  1196 1073 high
#7     0  1492  549 high
table(merged_data$risk)
#high  low 
#257  320

for(i in 1:577){
  merged_data$sample_name[i]<-paste(merged_data$sample_name[i],i,sep = "_")
}
row.names(merged_data)<-merged_data$sample_name
merged_data<-merged_data[,-1]
merged_data<-merged_data[,-59428]
dim(merged_data)
#[1]   577 59427
# the first 320 rows indicate the low-risk groups 
# the following 257 rows indicate the high-risk groups
dim(merged_data)
#[1]   577 59427
merged_data<-t(merged_data)
write.csv(merged_data, "Count_highlow_risk.csv", quote = F)

Group <- factor(rep(c("high","low"),times=c(257,320)),levels = c("low","high"))
table(Group)
#Group
#low high 
#320  257 
design <- model.matrix(~0+Group)
colnames(design)= levels(Group)
rownames(design)=colnames(merged_data)
write.csv(merged_data,"dataexpr.csv",quote = F)
#merged_data<-read.csv("dataexpr.csv",header=T,row.names=1)


library("limma")
library("tidyverse")
library("stringr")
library("edgeR")
table(is.na(dataexpr))
#FALSE 
#26682723
dge <- DGEList(counts=merged_data)
dge <- calcNormFactors(dge)

v <- voom(dge,design, normalize="quantile")
fit <- lmFit(v, design)

constrasts = paste(rev(levels(Group)),collapse = "-")
cont.matrix <- makeContrasts(contrasts=constrasts,levels = design) 
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)
DEG = topTable(fit2, coef=constrasts, n=Inf)
DEG = na.omit(DEG)
logFC_t=1
P.Value_t = 0.05
k1 = (DEG$P.Value < P.Value_t)&(DEG$logFC < -logFC_t)
k2 = (DEG$P.Value < P.Value_t)&(DEG$logFC > logFC_t)
change = ifelse(k1,"DOWN",ifelse(k2,"UP","stable"))
DEG$change <- change
DEG_limma <- DEG
table(DEG_limma$change)
#DOWN stable     UP 
#449  58776    202
write.csv(DEG_limma,"deg_High_LowRisk.csv",quote=F)

up_genes<-subset(DEG_limma, change == "UP")
up_list<-rownames(up_genes)
write.table(up_list, "up_genes.txt", sep = "\t", row.names = F, col.names = F, quote = F)

down_genes<-subset(DEG_limma, change == "DOWN")
down_list<-rownames(down_genes)
write.table(down_list, "down_genes.txt", sep = "\t", row.names = F, col.names = F, quote = F)

library("ggplot2")
deg<-DEG_limma
head(deg)
deg$color<-ifelse(deg$P.Value<0.05 & abs(deg$logFC)>logFC_t, ifelse(deg$logFC< -logFC_t, "blue", "red"),"gray")
color<-c(red = "red", gray = "gray", blue = "blue")
p<-ggplot(deg, aes(logFC, -log10(P.Value), col = color)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = color) +
  labs(x="log2 (fold change)", y = "-log10(P-Value)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
pdf("DEG_High_Low_VP.pdf")
p
dev.off()


#https://metascape.org/gp/index.html#/main/step1
#Figure 1. Bar graph of enriched terms across input gene lists, colored by p-values.


