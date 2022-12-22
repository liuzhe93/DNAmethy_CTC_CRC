setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/05_ImmuneLandscape/04_scRNA/GSE166555")
rm(list = ls()) 

library("ggplot2")
dat<-read.table("CRC_GSE166555_CellMetainfo_table.tsv", sep = "\t", header = T)
#############################malignancy##################################################
pdf("Mali_immu_stro.pdf")
ggplot(data = dat, mapping = aes(x = UMAP_1, y = UMAP_2, colour= Celltype..malignancy.)) + 
  geom_point(size = 0.5) + theme(legend.position = "top") + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))

dev.off()

#############################cell type####################################################
pdf("Celltype.pdf")
ggplot(data = dat, mapping = aes(x = UMAP_1, y = UMAP_2, colour= Celltype..major.lineage.)) + 
  geom_point(size = 0.5) + theme(legend.position = "top") + 
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))
dev.off()

############################Risk score##################################################
library("Seurat")
library("hdf5r")
library("tidyverse")
library("GSVA")
library("msigdbr")

dat<-Read10X_h5("CRC_GSE166555_expression.h5")
dim(dat)
#[1] 21753 66050
dat <- CreateSeuratObject(dat, project = "CRC", min.cells = 3)
dim(dat)
#[1] 21753 66050
dat <- NormalizeData(dat) %>% FindVariableFeatures() %>% ScaleData()

gene_set<-c("TRIP10", "NGFR", "SLC48A1", "SRMS")
gene_set<-as.data.frame(gene_set)
gene_set$y<-rep("Risk Factor",4)
colnames(gene_set)<-c("Metagene","Item")
list<- split(as.matrix(gene_set)[,1], gene_set[,2])

## 测试counts数据
expr <- as.matrix(dat@assays$RNA@counts)
system.time({res.counts = gsva(expr, list, method="ssgsea", parallel.sz=10)}) 
#   用户    系统    流逝 
#1615.07   12.11 1629.90

riskfactor_score<-t(res.counts)
RFS<-data.frame(riskfactor_score)
RFS$Cell<-row.names(RFS)
celltype_data<-read.table("CRC_GSE166555_CellMetainfo_table.tsv", sep = "\t", header = T)
dim(RFS)
#[1] 66050     2
dim(celltype_data)
#[1] 66050    10
input_data<-merge(celltype_data,RFS,by="Cell")
dim(input_data)
#[1] 66050    11
input_data$RiskFactorScore<-input_data$Risk.Factor
library("ggplot2")
pdf("RiskFactors_UMAP_score.pdf")
ggplot(data = input_data, mapping = aes(x = UMAP_1, y = UMAP_2, colour= RiskFactorScore)) + 
  geom_point(size = 0.5) + theme(legend.position = "top") +
  scale_color_gradient(low = "white", high = "red") +
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))
dev.off()
write.csv(input_data,"RiskFactors_data_scRNA.csv",quote=F,row.names = F)

#####################################distribution##########################################
RFS_data<-read.csv("RiskFactors_data_scRNA.csv")

library("ggplot2")
library("ggpubr")
library("ggstatsplot")
library("ggsci")

pdf("RiskFactors_score.pdf",height=4,width=10)
ggplot(data = RFS_data, mapping = aes(x = Celltype..major.lineage., y = RiskFactorScore, fill= Celltype..major.lineage.)) +
  geom_violin() + geom_boxplot(width=0.2)+
  #  scale_fill_jco()+geom_jitter(shape=16,size=2,position=position_jitter(0.2)) +
  stat_compare_means()+guides(fill=FALSE)+theme_classic()+
  geom_point(size = 1) + theme(legend.position = "top") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust =1, vjust = 1)) +
  scale_color_gradient(low = "white", high = "red")
dev.off()

#################################bubble plot#####################################
GLMS_data<-read.csv("RiskFactors_data_scRNA.csv")
table(GLMS_data$Celltype..major.lineage.)
#   B       CD4Tconv           CD8T             DC    Endothelial     Epithelial    Fibroblasts      Malignant           Mast 
#3400          11130           2809            718            649          17185           1636          12239            544 
#Mono/Macro Myofibroblasts         Plasma        Tprolif 
#2661            722          12011            346
GLMS_data$group<-GLMS_data$Celltype..major.lineage.

GLMS_B<-subset(GLMS_data,group=="B")
GLMS_B_avg<-mean(GLMS_B$RiskFactorScore)
GLMS_B_bg<-subset(GLMS_data,group!="B")
GLMS_B_avg_rest<-mean(GLMS_B_bg$RiskFactorScore)
log2_B<-log2(GLMS_B_avg/GLMS_B_avg_rest)

GLMS_CD4Tconv<-subset(GLMS_data,group=="CD4Tconv")
GLMS_CD4Tconv_avg<-mean(GLMS_CD4Tconv$RiskFactorScore)
GLMS_CD4Tconv_CD4Tconvg<-subset(GLMS_data,group!="CD4Tconv")
GLMS_CD4Tconv_avg_rest<-mean(GLMS_CD4Tconv_CD4Tconvg$RiskFactorScore)
log2_CD4Tconv<-log2(GLMS_CD4Tconv_avg/GLMS_CD4Tconv_avg_rest)

GLMS_CD8T<-subset(GLMS_data,group=="CD8T")
GLMS_CD8T_avg<-mean(GLMS_CD8T$RiskFactorScore)
GLMS_CD8T_bg<-subset(GLMS_data,group!="CD8T")
GLMS_CD8T_avg_rest<-mean(GLMS_CD8T_bg$RiskFactorScore)
log2_CD8T<-log2(GLMS_CD8T_avg/GLMS_CD8T_avg_rest)

GLMS_DC<-subset(GLMS_data,group=="DC")
GLMS_DC_avg<-mean(GLMS_DC$RiskFactorScore)
GLMS_DC_bg<-subset(GLMS_data,group!="DC")
GLMS_DC_avg_rest<-mean(GLMS_DC_bg$RiskFactorScore)
log2_DC<-log2(GLMS_DC_avg/GLMS_DC_avg_rest)

GLMS_Endothelial<-subset(GLMS_data,group=="Endothelial")
GLMS_Endothelial_avg<-mean(GLMS_Endothelial$RiskFactorScore)
GLMS_Endothelial_bg<-subset(GLMS_data,group!="Endothelial")
GLMS_Endothelial_avg_rest<-mean(GLMS_Endothelial_bg$RiskFactorScore)
log2_Endothelial<-log2(GLMS_Endothelial_avg/GLMS_Endothelial_avg_rest)

GLMS_Epithelial<-subset(GLMS_data,group=="Epithelial")
GLMS_Epithelial_avg<-mean(GLMS_Epithelial$RiskFactorScore)
GLMS_Epithelial_bg<-subset(GLMS_data,group!="Epithelial")
GLMS_Epithelial_avg_rest<-mean(GLMS_Epithelial_bg$RiskFactorScore)
log2_Epithelial<-log2(GLMS_Epithelial_avg/GLMS_Epithelial_avg_rest)

GLMS_Fibroblasts<-subset(GLMS_data,group=="Fibroblasts")
GLMS_Fibroblasts_avg<-mean(GLMS_Fibroblasts$RiskFactorScore)
GLMS_Fibroblasts_bg<-subset(GLMS_data,group!="Fibroblasts")
GLMS_Fibroblasts_avg_rest<-mean(GLMS_Fibroblasts_bg$RiskFactorScore)
log2_Fibroblasts<-log2(GLMS_Fibroblasts_avg/GLMS_Fibroblasts_avg_rest)

GLMS_Malignant<-subset(GLMS_data,group=="Malignant")
GLMS_Malignant_avg<-mean(GLMS_Malignant$RiskFactorScore)
GLMS_Malignant_bg<-subset(GLMS_data,group!="Malignant")
GLMS_Malignant_avg_rest<-mean(GLMS_Malignant_bg$RiskFactorScore)
log2_Malignant<-log2(GLMS_Malignant_avg/GLMS_Malignant_avg_rest)

GLMS_Mast<-subset(GLMS_data,group=="Mast")
GLMS_Mast_avg<-mean(GLMS_Mast$RiskFactorScore)
GLMS_Mast_bg<-subset(GLMS_data,group!="Mast")
GLMS_Mast_avg_rest<-mean(GLMS_Mast_bg$RiskFactorScore)
log2_Mast<-log2(GLMS_Mast_avg/GLMS_Mast_avg_rest)

GLMS_MonoMacro<-subset(GLMS_data,group=="Mono/Macro")
GLMS_MonoMacro_avg<-mean(GLMS_MonoMacro$RiskFactorScore)
GLMS_MonoMacro_bg<-subset(GLMS_data,group!="Mono/Macro")
GLMS_MonoMacro_avg_rest<-mean(GLMS_MonoMacro_bg$RiskFactorScore)
log2_MonoMacro<-log2(GLMS_MonoMacro_avg/GLMS_MonoMacro_avg_rest)

GLMS_Myofibroblasts<-subset(GLMS_data,group=="Myofibroblasts")
GLMS_Myofibroblasts_avg<-mean(GLMS_Myofibroblasts$RiskFactorScore)
GLMS_Myofibroblasts_bg<-subset(GLMS_data,group!="Myofibroblasts")
GLMS_Myofibroblasts_avg_rest<-mean(GLMS_Myofibroblasts_bg$RiskFactorScore)
log2_Myofibroblasts<-log2(GLMS_Myofibroblasts_avg/GLMS_Myofibroblasts_avg_rest)

GLMS_Plasma<-subset(GLMS_data,group=="Plasma")
GLMS_Plasma_avg<-mean(GLMS_Plasma$RiskFactorScore)
GLMS_Plasma_bg<-subset(GLMS_data,group!="Plasma")
GLMS_Plasma_avg_rest<-mean(GLMS_Plasma_bg$RiskFactorScore)
log2_Plasma<-log2(GLMS_Plasma_avg/GLMS_Plasma_avg_rest)

GLMS_Tprolif<-subset(GLMS_data,group=="Tprolif")
GLMS_Tprolif_avg<-mean(GLMS_Tprolif$RiskFactorScore)
GLMS_Tprolif_bg<-subset(GLMS_data,group!="Tprolif")
GLMS_Tprolif_avg_rest<-mean(GLMS_Tprolif_bg$RiskFactorScore)
log2_Tprolif<-log2(GLMS_Tprolif_avg/GLMS_Tprolif_avg_rest)

data.matrix<-matrix(c(log2_B,"B",log2_CD4Tconv,"CD4Tconv",log2_CD8T,"CD8T",log2_DC,"DC",log2_Endothelial,"Endothelial",
                      log2_Epithelial,"Epithelial",log2_Fibroblasts,"Fibroblasts",log2_Malignant,"Malignant",
                      log2_Mast,"Mast",log2_MonoMacro,"Mono/Macro",log2_Myofibroblasts,"Myofibroblasts",log2_Plasma,"Plasma",
                      log2_Tprolif,"Tprolif"),c(13,2),byrow = T)
data<-data.frame(data.matrix)
colnames(data)<-c("log2FC","CellType")
data$log2FC<-as.numeric(data$log2FC)
library("ggplot2")

pdf("RiskFactor_enrich.pdf",width=4,height=4)
ggplot(data,aes(log2FC,CellType))+geom_point(aes(color=log2FC),size=6) +
  scale_color_gradient2(low = "blue", mid = "white", high =  "red") +
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))
dev.off()



