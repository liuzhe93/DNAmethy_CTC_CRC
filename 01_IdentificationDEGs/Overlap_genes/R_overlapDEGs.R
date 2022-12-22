setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/Overlap_genes")
rm(list=ls())

deg_CTC<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/CTC_PT/deg_CTC_PT.csv",header=T)
deg_metastic<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/Metastatic_nonMetastatic/deg_Me_nonMe.csv",header=T)

A<-subset(deg_CTC,g!="stable")
B<-subset(deg_metastic,change!="stable")
library(VennDiagram)
Length_A<-length(A$X)
Length_B<-length(B$X)
Length_AB<-length(intersect(A$X,B$X))
Length_A
#[1] 9363
Length_B
#[1] 11702
Length_AB
#[1] 1573
T<-venn.diagram(list("CTC_PT"=A$X,"Me_NonMe"=B$X),filename=NULL,lwd=1,lty=2, ,col=c('red','green'),fill=c('red','green'),cat.col=c('red','green'))
pdf("OverlapDEGs.pdf")
grid.draw(T)
dev.off()

######################################################CTC_PT####################################################################################
merged_AB<-merge(A, B, by = "X")
merged_AB_sorted<-merged_AB[order(merged_AB$g, decreasing = TRUE),]
library(pheatmap)
genes<-merged_AB_sorted$X
FPKM_CTC_PT<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/CTC_PT/FPKM_CTC_PT.csv",row.names = 1, header = TRUE)
colnames(FPKM_CTC_PT) = c(paste("PT", 1:9, sep = ""),paste("CTC", 1:18, sep = ""))
FPKM_CTC_PT_selected<-FPKM_CTC_PT[genes,]
dim(FPKM_CTC_PT_selected)
#[1] 1573  27
annotation_col = data.frame( CellType = factor(c(rep("PT",9), rep("CTC",18)),c("PT","CTC")))
rownames(annotation_col) = c(paste("PT", 1:9, sep = ""),paste("CTC", 1:18, sep = ""))
annotation_row = data.frame( GeneClass = factor(c(rep("UP",1259), rep("DOWN",314)),c("UP","DOWN")))
rownames(annotation_row) = merged_AB_sorted$X
ann_colors = list(CellType = c(CTC = "red", PT = "blue"), GeneClass = c(UP = "red", DOWN = "blue") )
pdf("heatmap_CTC_PT.pdf")
pheatmap(FPKM_CTC_PT_selected,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F, 
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col, annotation_row = annotation_row,
         angle_col = "45", annotation_colors = ann_colors, main = "Title", gaps_row = 1259, gaps_col = 9, fontsize_row = 4)
dev.off()
############################################Metastatic_Non-metastatic###############################################################
merged_AB_sorted<-merged_AB[order(merged_AB$change, decreasing = TRUE),]
library(pheatmap)
genes<-merged_AB_sorted$X
FPKM_Me_noME<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/Data_Collection/TCGA-COAD/mRNA_exp/FPKM.csv")
library(rtracklayer)
x = rtracklayer::import("gencode.v36.annotation.gtf.gz")
x2 = as.data.frame(x)
tj = as.data.frame(table(x2$type))
tj
#            Var1    Freq
#1           gene   60660
#2     transcript  232117
#3           exon 1429877
#4            CDS  790200
#5    start_codon   90112
#6     stop_codon   82924
#7            UTR  328848
#8 Selenocysteine     117
anno<-subset(x2,type == "gene", select = c("gene_name","gene_id"))
FPKM_anno<-merge(anno,FPKM_Me_noME,by.x="gene_id",by.y="X")
FPKM_anno<-FPKM_anno[,-1]
FPKM_anno<-aggregate(FPKM_anno,by=list(FPKM_anno$gene_name),FUN=median)
#FPKM_anno<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/01_univariate_Cox_regression_analysis/FPKM_anno_uniq.csv",
#                    row.names=1)
#rownames(FPKM_anno)<-FPKM_anno$Group.1
#FPKM_anno<-FPKM_anno[,c(-1,-2)]
dim(FPKM_anno)
#[1] 59427   521
FPKM_anno_t<-t(FPKM_anno)
FPKM_anno_t<-as.data.frame(FPKM_anno_t)
FPKM_anno_t$sample_name<-gsub("\\.","-",rownames(FPKM_anno_t))
FPKM_anno_t$sample_name<-substr(FPKM_anno_t$sample_name,1,12)
FPKM_anno_t$full_name<-rownames(FPKM_anno_t)
result_time<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/Metastatic_nonMetastatic/result_time.csv",header=TRUE)
result_time$sample_type<-ifelse(result_time$ajcc_pathologic_m=="M0","NonMetastatic",
                                ifelse(result_time$ajcc_pathologic_m=="M1","Metastatic",
                                       ifelse(result_time$ajcc_pathologic_m=="M1a","Metastatic",
                                              ifelse(result_time$ajcc_pathologic_m=="M1b","Metastatic",
                                                     "NotSure"))))
FPKM_selected<-merge(FPKM_anno_t,result_time,by.x="sample_name",by.y="submitter_id")
FPKM_nonme<-subset(FPKM_selected,sample_type=="NonMetastatic")
FPKM_nonme<-FPKM_nonme[,1:59429]
dim(FPKM_nonme)
#[1]   376 59429
# there are 376 non-metastatic samples
FPKM_me<-subset(FPKM_selected,sample_type=="Metastatic")
FPKM_me<-FPKM_me[,1:59429]
dim(FPKM_me)
#[1]    73 59429
# there are 73 metastatic samples
dat<-rbind(FPKM_nonme,FPKM_me)
dat<-dat[,c(59429,1:59428)]
dat<-dat[,-2]
dat_t<-t(dat)
colnames(dat_t)<-dat_t[1,]
dat_t<-dat_t[-1,]
dat_t<-as.data.frame(dat_t)
dat_t$gene_name<-row.names(dat_t)
data_numeric <- apply(dat_t[,1:(dim(dat_t)[2]-1)],2,as.numeric)
dataexpr <- data.frame(data_numeric, dat_t[,dim(dat_t)[2]])
row.names(dataexpr)<-dataexpr[,dim(dat_t)[2]]
dataexpr<-dataexpr[,-dim(dat_t)[2]]
dim(dataexpr)
write.csv(dataexpr,"FPKM_mRNAExp.Me_NoMe_selected.csv",quote=F)
colnames(dataexpr) = c(paste("non_metastatic", 1:376, sep = ""),paste("metastatic", 1:73, sep = ""))
FPKM_Me_NoMe_selected<-dataexpr[genes,]
dim(FPKM_Me_NoMe_selected)
#[1] 1573  449
annotation_col = data.frame( CellType = factor(c(rep("non_metastatic",376), rep("metastatic",73)),c("non_metastatic","metastatic")))
rownames(annotation_col) = c(paste("non_metastatic", 1:376, sep = ""),paste("metastatic", 1:73, sep = ""))
annotation_row = data.frame( GeneClass = factor(c(rep("UP",801), rep("DOWN",772)),c("UP","DOWN")))
rownames(annotation_row) = merged_AB_sorted$X
ann_colors = list(CellType = c(metastatic = "red", non_metastatic = "blue"), GeneClass = c(UP = "red", DOWN = "blue") )
pdf("heatmap_Me_NoMe.pdf")
pheatmap(FPKM_Me_NoMe_selected,scale = "row", clustering_distance_rows = "correlation", cluster_rows = F, cluster_col = F, 
         color = colorRampPalette(c("blue", "white", "red"))(50), annotation_col = annotation_col, annotation_row = annotation_row,
         angle_col = "45", annotation_colors = ann_colors, main = "Title", gaps_row = 801, gaps_col = 376, fontsize_row = 4, fontsize_col = 4)
dev.off()
####################################################################################################################################
methy<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/Data_Collection/TCGA-COAD/methylation/TCGA_coad.met_uniq.csv",header = TRUE, row.names = 1)
dataexpr<-read.csv("FPKM_mRNAExp.Me_NoMe_selected.csv",header=TRUE,row.names = 1)
merged_AB<-merge(A, B, by = "X")
dim(merged_AB)
#[1] 1573  15
dim(methy)
#[1] 24311   351
genes<-merged_AB$X
methy_selected<-methy[genes,]
methy_selected_rmNA<-na.omit(methy_selected)
FPKM_selected<-dataexpr
FPKM_selected<-FPKM_selected[genes,]
dim(FPKM_selected)
#[1] 1573 449
dim(methy_selected_rmNA)
#[1] 1217 351
colnames(FPKM_selected)<-substr(colnames(FPKM_selected),1,16)
colnames(methy_selected_rmNA)<-substr(colnames(methy_selected_rmNA),1,16)
sameSamples<-intersect(colnames(FPKM_selected),colnames(methy_selected_rmNA))
length(sameSamples)
#[1] 261
FPKM_selected1=FPKM_selected[,sameSamples]
methy_selected1=methy_selected_rmNA[,sameSamples]
sameGenes<-intersect(rownames(FPKM_selected),rownames(methy_selected_rmNA))
length(sameGenes)
#[1] 1146
FPKM_selected2=FPKM_selected1[sameGenes,]
methy_selected2=methy_selected1[sameGenes,]

corFilter=0             #相关系数过滤标准
pvalueFilter=0.05         #p值过滤标准
outTab=data.frame()
for(i in row.names(FPKM_selected2)){
  if(sd(FPKM_selected2[i,])>1){
#    for(j in row.names(methy_selected2)){
      x=as.numeric(FPKM_selected2[i,])
      y=as.numeric(methy_selected2[i,])
      corT=cor.test(x,y)
      cor=corT$estimate
      pvalue=corT$p.value
#      if((cor>corFilter) & (pvalue<pvalueFilter)){
#        outTab=rbind(outTab,cbind(FPKM=i,methy=j,cor,pvalue,Regulation="postive"))
#      }
      if((cor< -corFilter) & (pvalue<pvalueFilter)){
        outTab=rbind(outTab,cbind(FPKM=i,methy=i,cor,pvalue,Regulation="negative"))
#      }
    }
  }
}
write.csv(outTab,file="corResult.csv",quote=F,row.names=F)        #输出相关性结果
dim(outTab)
#[1] 166  5

outTab_order<-outTab[order(outTab$cor,decreasing = TRUE),]
top6<-head(outTab_order,n=6)

library("ggExtra")
library("ggplot2")
library("ggpubr")
library("ggpmisc")
theme_set(ggpubr::theme_pubr()+
            theme(legend.position = "top"))
#####################################################  ELF5 ##############################################################
xval<-t(methy_selected2["ELF5",])
yval<-t(FPKM_selected2["ELF5",])
data_ELF5<-cbind(xval,yval)
colnames(data_ELF5)<-c("Methylation","Expression")
data_ELF5<-as.data.frame(data_ELF5)
pdf("ELF5.pdf")
p<-ggscatter(data_ELF5,x="Methylation",y="Expression",
          add = "reg.line", conf.int = T,
          add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
#####################################################  PIGR ##############################################################
xval<-t(methy_selected2["PIGR",])
yval<-t(FPKM_selected2["PIGR",])
data_PIGR<-cbind(xval,yval)
colnames(data_PIGR)<-c("Methylation","Expression")
data_PIGR<-as.data.frame(data_PIGR)
pdf("PIGR.pdf")
p<-ggscatter(data_PIGR,x="Methylation",y="Expression",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
#####################################################  CLDN2 ##############################################################
xval<-t(methy_selected2["CLDN2",])
yval<-t(FPKM_selected2["CLDN2",])
data_CLDN2<-cbind(xval,yval)
colnames(data_CLDN2)<-c("Methylation","Expression")
data_CLDN2<-as.data.frame(data_CLDN2)
pdf("CLDN2.pdf")
p<-ggscatter(data_CLDN2,x="Methylation",y="Expression",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
#####################################################  PLA2G3 ##############################################################
xval<-t(methy_selected2["PLA2G3",])
yval<-t(FPKM_selected2["PLA2G3",])
data_PLA2G3<-cbind(xval,yval)
colnames(data_PLA2G3)<-c("Methylation","Expression")
data_PLA2G3<-as.data.frame(data_PLA2G3)
pdf("PLA2G3.pdf")
p<-ggscatter(data_PLA2G3,x="Methylation",y="Expression",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
#####################################################  ARHGAP30 ##############################################################
xval<-t(methy_selected2["ARHGAP30",])
yval<-t(FPKM_selected2["ARHGAP30",])
data_ARHGAP30<-cbind(xval,yval)
colnames(data_ARHGAP30)<-c("Methylation","Expression")
data_ARHGAP30<-as.data.frame(data_ARHGAP30)
pdf("ARHGAP30.pdf")
p<-ggscatter(data_ARHGAP30,x="Methylation",y="Expression",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
#####################################################  KRT23 ##############################################################
xval<-t(methy_selected2["KRT23",])
yval<-t(FPKM_selected2["KRT23",])
data_KRT23<-cbind(xval,yval)
colnames(data_KRT23)<-c("Methylation","Expression")
data_KRT23<-as.data.frame(data_KRT23)
pdf("KRT23.pdf")
p<-ggscatter(data_KRT23,x="Methylation",y="Expression",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
#####################################################  SLAMF8 ##############################################################
xval<-t(methy_selected2["SLAMF8",])
yval<-t(FPKM_selected2["SLAMF8",])
data_SLAMF8<-cbind(xval,yval)
colnames(data_SLAMF8)<-c("Methylation","Expression")
data_SLAMF8<-as.data.frame(data_SLAMF8)
pdf("SLAMF8.pdf")
p<-ggscatter(data_SLAMF8,x="Methylation",y="Expression",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
#####################################################  CAPS ##############################################################
xval<-t(methy_selected2["CAPS",])
yval<-t(FPKM_selected2["CAPS",])
data_CAPS<-cbind(xval,yval)
colnames(data_CAPS)<-c("Methylation","Expression")
data_CAPS<-as.data.frame(data_CAPS)
pdf("CAPS.pdf")
p<-ggscatter(data_CAPS,x="Methylation",y="Expression",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
#####################################################  AZGP1 ##############################################################
xval<-t(methy_selected2["AZGP1",])
yval<-t(FPKM_selected2["AZGP1",])
data_AZGP1<-cbind(xval,yval)
colnames(data_AZGP1)<-c("Methylation","Expression")
data_AZGP1<-as.data.frame(data_AZGP1)
pdf("AZGP1.pdf")
p<-ggscatter(data_AZGP1,x="Methylation",y="Expression",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
#####################################################  UBE2L6 ##############################################################
xval<-t(methy_selected2["UBE2L6",])
yval<-t(FPKM_selected2["UBE2L6",])
data_UBE2L6<-cbind(xval,yval)
colnames(data_UBE2L6)<-c("Methylation","Expression")
data_UBE2L6<-as.data.frame(data_UBE2L6)
pdf("UBE2L6.pdf")
p<-ggscatter(data_UBE2L6,x="Methylation",y="Expression",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
#####################################################  LYZ ##############################################################
xval<-t(methy_selected2["LYZ",])
yval<-t(FPKM_selected2["LYZ",])
data_LYZ<-cbind(xval,yval)
colnames(data_LYZ)<-c("Methylation","Expression")
data_LYZ<-as.data.frame(data_LYZ)
pdf("LYZ.pdf")
p<-ggscatter(data_LYZ,x="Methylation",y="Expression",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
#####################################################  POLR1D ##############################################################
xval<-t(methy_selected2["POLR1D",])
yval<-t(FPKM_selected2["POLR1D",])
data_POLR1D<-cbind(xval,yval)
colnames(data_POLR1D)<-c("Methylation","Expression")
data_POLR1D<-as.data.frame(data_POLR1D)
pdf("POLR1D.pdf")
p<-ggscatter(data_POLR1D,x="Methylation",y="Expression",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
#############################################################################################################################
head(outTab_order,n=12)
#FPKM   methy                cor               pvalue Regulation
#cor35      ELF5     ELF5 -0.630783909024941 2.29447801189009e-30   negative
#cor102     PIGR     PIGR -0.620875415110866 3.26657716382768e-29   negative
#cor20     CLDN2    CLDN2 -0.574943415502891 2.29934475613349e-24   negative
#cor105   PLA2G3   PLA2G3 -0.560302781652734 5.64730482190669e-23   negative
#cor1   ARHGAP30 ARHGAP30 -0.558252854253986  8.7307538722766e-23   negative
#cor66     KRT23    KRT23 -0.542492280773521 2.25449332229405e-21   negative
#cor128   SLAMF8   SLAMF8 -0.508584624660344 1.41659936935576e-18   negative
#cor8       CAPS     CAPS -0.484011466742049 9.83438047971729e-17   negative
#cor3      AZGP1    AZGP1 -0.481718573376605  1.4359739280425e-16   negative
#cor159   UBE2L6   UBE2L6 -0.466017176021206 1.78014570189869e-15   negative
#cor72       LYZ      LYZ  -0.46253922556725  3.0559063375675e-15   negative
#cor106   POLR1D   POLR1D -0.447845273861299 2.80373034592671e-14   negative
# In total, there are 166 genes were kept.
# In this sutdy, we only selected six genes for visualization, ELF5, CLDN2, GNG4, CD3D, TNNC2, and PLA2G2A.
# The reason is that the distribution of their data distribution.

