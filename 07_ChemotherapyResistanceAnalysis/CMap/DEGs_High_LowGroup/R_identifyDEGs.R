setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/07_ChemotherapyResistanceAnalysis/CMap/DEGs_High_LowGroup")
rm(list = ls())

risk<-read.table("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/04_multivariate_Cox_regression_analysis/risk.txt",
                 header=T,row.names = 1)
risk$sample_name<-row.names(risk)
risk$sample_name<-substr(risk$sample_name,1,12)
dim(risk)
#[1] 521   9
head(risk)

counts_TCGA <- read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/Metastatic_nonMetastatic/dataexpr.csv",
                                header = T,check.names = FALSE)
# 定义行名
rownames(counts_TCGA) <- counts_TCGA[,1]
counts_TCGA <- counts_TCGA[,-1]
# 先看一下矩阵长啥样，心里有个数：每一行是一个基因，每一列是一个样本
#View(counts_TCGA)
dim(counts_TCGA)
# [1] 59427   449
counts_TCGA_t<-t(counts_TCGA)
counts_TCGA_t<-as.data.frame(counts_TCGA_t)
counts_TCGA_t$sample_name<-row.names(counts_TCGA_t)
counts_TCGA_t$sample_name<-substr(counts_TCGA_t$sample_name, 1, 12)
counts_TCGA_t$sample_name<-gsub("\\.", "-", counts_TCGA_t$sample_name)

merged_data<-merge(counts_TCGA_t, risk, by = "sample_name")
merged_data<-subset(merged_data, select = c(1:59428,59436))
sample_id2type<-cbind(merged_data$sample_name, merged_data$risk)
colnames(sample_id2type)<-c("sample_name", "risk_type")
merged_data_high<-subset(merged_data, risk == "high")
merged_data_low<-subset(merged_data, risk == "low")
dim(merged_data_high)
#[1]   257 59429
dim(merged_data_low)
#[1]   320 59429
merged_data_high<-merged_data_high[,-59429]
merged_data_low<-merged_data_low[,-59429]

mydata<-rbind(merged_data_high, merged_data_low)
dim(mydata)
#[1]   577 59428
for(i in 1:nrow(mydata)){
  mydata$sample_name[i] <- paste(mydata$sample_name[i], i, sep = "_")
}
row.names(mydata)<-mydata$sample_name
mydata<-mydata[,-1]
coad_counts<-t(mydata)
dim(coad_counts)
#[1] 59427   577
coad_counts<-as.data.frame(coad_counts)
Group <- factor(rep(c("high","low"),times=c(257,320)),levels = c("low","high"))
table(Group)
#Group
#low high 
#320  257
design <- model.matrix(~0+Group)
colnames(design)= levels(Group)
rownames(design)=colnames(coad_counts)

library("limma")
library("tidyverse")
library("stringr")
library("edgeR")
table(is.na(coad_counts))
#FALSE 
#34289379
dge <- DGEList(counts=coad_counts)
dge <- calcNormFactors(dge)
v <- voom(dge,design, normalize="quantile")
fit <- lmFit(v, design)

constrasts = paste(rev(levels(Group)),collapse = "-")
cont.matrix <- makeContrasts(contrasts=constrasts,levels = design) 
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)
DEG = topTable(fit2, coef=constrasts, n=Inf)
DEG = na.omit(DEG)
#logFC_t=0.5849625
logFC_t=1
P.Value_t = 0.05
k1 = (DEG$P.Value < P.Value_t)&(DEG$logFC < -logFC_t)
k2 = (DEG$P.Value < P.Value_t)&(DEG$logFC > logFC_t)
change = ifelse(k1,"DOWN",ifelse(k2,"UP","stable"))
DEG$change <- change
DEG_limma <- DEG
table(DEG_limma$change)
#DOWN stable     UP 
# 449  58776    202 
write.csv(DEG_limma,"deg_High_Low.csv",quote=F)

#######clue query: https://clue.io/query###############
#ref: https://zhuanlan.zhihu.com/p/520458579
#输入的差异基因（上调和下调基因10-150个），只有蓝色打钩的基因才是可用的，我们可以把红色和空白圈圈的基因去除。最后点submit。



