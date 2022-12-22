setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/06_GenomicVariation/CNV")
rm(list = ls()) 

#ref: https://www.jianshu.com/p/b8351a86a40b


#####################################下载并整理TCGA的CNV数据##########################################
library("TCGAbiolinks")
cancer_type <- "TCGA-COAD"
clinical <- GDCquery_clinic(project = cancer_type, type = "clinical")
library("SummarizedExperiment")
data_type <- "Masked Copy Number Segment"
data_category <- "Copy Number Variation"
workflow_type <- "DNAcopy"
query_coad_cnv <- GDCquery(project = cancer_type, 
                           data.category = data_category,
                           data.type = data_type,
                           workflow.type = workflow_type)
#GDCdownload(query_coad_cnv, method = "api", files.per.chunk = 100)
#CNV_files <- GDCprepare(query = query_coad_cnv, save = TRUE, save.filename = "CNV_TCGA_COAD.rda")

cnv <- load("CNV_TCGA_COAD.rda")
coad_seg <- eval(parse(text = cnv))
coad_seg <- coad_seg[,-1]
coad_seg <- coad_seg[,c('Sample','Chromosome','Start','End','Num_Probes','Segment_Mean')]
#tumor_seg <- coad_seg[substr(coad_seg$Sample,14,15) == "01",]
#write.table(tumor_seg, file = "COAD_CNV.txt", sep = "\t", quote = F, row.names = F)
write.table(coad_seg, file = "COAD_CNV.txt", sep = "\t", quote = F, row.names = F)

#################################Marker文件的下载以及格式格式调整#####################################
hg_marker_file <- read.delim("snp6.na35.remap.hg38.subset.txt.gz")
View(hg_marker_file)
hg_marker_file <- hg_marker_file[hg_marker_file$freqcnv == "FALSE",]
hg_marker_file <- hg_marker_file[,c(1,2,3)]
write.table(hg_marker_file, "hg_marker_file.txt", sep = "\t", col.names = TRUE, row.names = F)

####################################运行GISTIC（online版）#############################################
#high risk score group 和 low risk score group 分别运行 分别得出结果

#############################Violin plot: amplification and deletion#################################################
setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/06_GenomicVariation_CNV/Results_GenePattern")
rm(list = ls()) 
library("maftools")
coad.gistic <- readGistic(gisticAllLesionsFile="all_lesions.conf_90.txt", gisticAmpGenesFile="amp_genes.conf_90.txt", 
                          gisticDelGenesFile="del_genes.conf_90.txt", gisticScoresFile="scores.gistic", isTCGA=TRUE)

gisticChromPlot(gistic = coad.gistic, ref.build = "hg38")


##############################################GenePattern#################################################################
setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/06_GenomicVariation_CNV/CNV")
rm(list = ls()) 
tumor_seg<-read.table("COAD_CNV.txt", sep = "\t", header = T)
tumor_seg$sample_name<-substr(tumor_seg$Sample,1,12)
risk<-read.table("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/04_multivariate_Cox_regression_analysis/risk.txt",
                 sep="\t", header=T, row.names = 1)
risk$sample_name<-row.names(risk)
risk$sample_name<-substr(risk$sample_name,1,12)
head(risk)
merged_data<-merge(tumor_seg, risk, by = "sample_name")
mydata_high<-subset(merged_data, risk == "high")
mydata_high<-subset(mydata_high, select = c( "Sample","Chromosome","Start","End","Num_Probes","Segment_Mean"))
mydata_high_dup<-mydata_high[!duplicated(mydata_high, fromLast=TRUE), ] 
write.table(mydata_high_dup, file = "COAD_CNV_high.txt", sep = "\t", quote = F, row.names = F)
mydata_low<-subset(merged_data, risk == "low")
mydata_low<-subset(mydata_low, select = c( "Sample","Chromosome","Start","End","Num_Probes","Segment_Mean"))
mydata_low_dup<-mydata_low[!duplicated(mydata_low, fromLast=TRUE), ] 
write.table(mydata_low_dup, file = "COAD_CNV_low.txt", sep = "\t", quote = F, row.names = F)

setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/06_GenomicVariation_CNV/Results_GenePattern/high_results/")
rm(list = ls()) 
library("maftools")
coad.gistic <- readGistic(gisticAllLesionsFile="all_lesions.conf_90.txt", gisticAmpGenesFile="amp_genes.conf_90.txt", 
                          gisticDelGenesFile="del_genes.conf_90.txt", gisticScoresFile="scores.gistic", isTCGA=TRUE)
mydata_high<-getSampleSummary(coad.gistic)
pdf("HighGroup_CNV_GenomePlot.pdf",width = 10, height = 6)
gisticChromPlot(gistic = coad.gistic, markBands = "all")
dev.off()
pdf("HighGroup_BubblePlot.pdf")
gisticBubblePlot(gistic = coad.gistic)
dev.off()

setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/06_GenomicVariation_CNV/Results_GenePattern/low_results/")
coad.gistic <- readGistic(gisticAllLesionsFile="all_lesions.conf_90.txt", gisticAmpGenesFile="amp_genes.conf_90.txt", 
                          gisticDelGenesFile="del_genes.conf_90.txt", gisticScoresFile="scores.gistic", isTCGA=TRUE)
mydata_low<-getSampleSummary(coad.gistic)
pdf("LowGroup_CNV_GenomePlot.pdf",width = 10, height = 6)
gisticChromPlot(gistic = coad.gistic, markBands = "all")
dev.off()
pdf("LowGroup_BubblePlot.pdf")
gisticBubblePlot(gistic = coad.gistic)
dev.off()

setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/06_GenomicVariation_CNV/Results_GenePattern/amp_del")
dim(mydata_high)
#[1] 253   4
dim(mydata_low)
#[1] 230   4
bind_mydata<-cbind(mydata_high, mydata_low)
bind_mydata<-bind_mydata[1:230,]
bind_mydata<-as.data.frame(bind_mydata)
colnames(bind_mydata)<-c("Barcode_high","Amp_high","Del_high","Total_high","Barcode_low","Amp_low","Del_low","Total_low")

###########################amplification boxplot##########################################################################
bind_mydata_amp<-subset(bind_mydata, select = c("Barcode_high","Amp_high","Amp_low"))
for(i in 1:nrow(bind_mydata_amp)){
  bind_mydata_amp$sample_name[i] <- paste(bind_mydata_amp$Barcode_high[i], i, sep = "_")
}
rownames(bind_mydata_amp)<-bind_mydata_amp$sample_name
bind_mydata_amp<-subset(bind_mydata_amp, select = c("Amp_high", "Amp_low"))
head(bind_mydata_amp)
#               Amp_high Amp_low
#TCGA-A6-2671_1     1164     710
#TCGA-A6-6650_2     2181    1000
#TCGA-A6-4107_3     1485    1238
#TCGA-A6-2681_4     1640     562
#TCGA-G4-6317_5     1269    1190
#TCGA-NH-A8F7_6     1468     418
dat <- bind_mydata_amp %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = CNV_Type,value = Mutation,-Sample)
dat$Mutation<-log2(dat$Mutation)
library(tinyarray)
library(tidyverse)
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
library(ggpubr)

pdf("AMP_boxplot.pdf")
p <- ggplot(dat, aes(x = CNV_Type, y = Mutation, fill = CNV_Type)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge()) +
  scale_fill_aaas() +
  theme_classic(base_size = 16)+
  labs(x = "", y = 'Amp log2 Mutation Counts') +
  stat_compare_means() +
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ stat_compare_means(aes(group = CNV_Type,label = ..p.signif..),method = "kruskal.test")
p
dev.off()
###########################deletion boxplot##########################################################################
bind_mydata_del<-subset(bind_mydata, select = c("Barcode_high","Del_high","Del_low"))
for(i in 1:nrow(bind_mydata_del)){
  bind_mydata_del$sample_name[i] <- paste(bind_mydata_del$Barcode_high[i], i, sep = "_")
}
rownames(bind_mydata_del)<-bind_mydata_del$sample_name
bind_mydata_del<-subset(bind_mydata_del, select = c("Del_high", "Del_low"))
head(bind_mydata_del)
#               Del_high Del_low
#TCGA-A6-2671_1     3942    9996
#TCGA-A6-6650_2     2571    9126
#TCGA-A6-4107_3     2942    6621
#TCGA-A6-2681_4     2622    6514
#TCGA-G4-6317_5     2550    5480
#TCGA-NH-A8F7_6     2239    5599
dat <- bind_mydata_del %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = CNV_Type,value = Mutation,-Sample)
dat$Mutation<-log2(dat$Mutation)
library(tinyarray)
library(tidyverse)
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
library(ggpubr)

pdf("DEL_boxplot.pdf")
p <- ggplot(dat, aes(x = CNV_Type, y = Mutation, fill = CNV_Type)) +
  geom_boxplot(position = position_dodge(0.8)) +
  geom_point(position = position_jitterdodge()) +
  scale_fill_aaas() +
  theme_classic(base_size = 16)+
  labs(x = "", y = 'Del log2 Mutation Counts') +
  stat_compare_means() +
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ stat_compare_means(aes(group = CNV_Type,label = ..p.signif..),method = "kruskal.test")
p
dev.off()


