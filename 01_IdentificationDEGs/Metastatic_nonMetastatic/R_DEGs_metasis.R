setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/Metastatic_nonMetastatic")
rm(list = ls())

surv_data<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/Data_Collection/TCGA-COAD/clinical/TCGAbiolinks-COAD-clinical.csv")
surv_data<-subset(surv_data,select = c("submitter_id","ajcc_pathologic_stage","ajcc_pathologic_t","ajcc_pathologic_n","ajcc_pathologic_m",
                                     "gender","vital_status","age_at_index","days_to_last_follow_up","days_to_death"))
surv_data<-surv_data
surv_data$os<-ifelse(surv_data$vital_status=='Alive',surv_data$days_to_last_follow_up,surv_data$days_to_death)
library(dplyr)
surv_selected<-surv_data %>%
  select(submitter_id,os,vital_status,gender,age_at_index,ajcc_pathologic_stage,ajcc_pathologic_t,ajcc_pathologic_m,ajcc_pathologic_n)
#0=alive, 1=dead.
surv_selected$fustat<-ifelse(surv_selected$vital_status=='Alive',0,1)
surv_selected$futime<-surv_selected$os/365
result_time<-surv_selected %>%
  select(submitter_id,futime,fustat,gender,age_at_index,ajcc_pathologic_stage,ajcc_pathologic_m)
write.csv(result_time,"result_time.csv",row.names=F)
table(result_time$ajcc_pathologic_m)

counts<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/Data_Collection/TCGA-COAD/mRNA_exp/count.csv")
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
counts_anno<-merge(anno,counts,by.x="gene_id",by.y="X")
counts_anno<-counts_anno[,-1]
counts_anno<-aggregate(counts_anno,by=list(counts_anno$gene_name),FUN=median)
rownames(counts_anno)<-counts_anno$Group.1
counts_anno<-counts_anno[,c(-1,-2)]
dim(counts_anno)
#[1] 59427   521
counts_anno_t<-t(counts_anno)
counts_anno_t<-as.data.frame(counts_anno_t)
counts_anno_t$sample_name<-gsub("\\.","-",rownames(counts_anno_t))
counts_anno_t$sample_name<-substr(counts_anno_t$sample_name,1,12)
result_time$sample_type<-ifelse(result_time$ajcc_pathologic_m=="M0","NonMetastatic",
                                ifelse(result_time$ajcc_pathologic_m=="M1","Metastatic",
                                       ifelse(result_time$ajcc_pathologic_m=="M1a","Metastatic",
                                              ifelse(result_time$ajcc_pathologic_m=="M1b","Metastatic",
                                       "NotSure"))))
counts_selected<-merge(counts_anno_t,result_time,by.x="sample_name",by.y="submitter_id")
counts_nonme<-subset(counts_selected,sample_type=="NonMetastatic")
counts_nonme<-counts_nonme[,1:59428]
counts_me<-subset(counts_selected,sample_type=="Metastatic")
counts_me<-counts_me[,1:59428]
dim(counts_me)
#[1]    73 59428
# there are 73 metastatic samples
dim(counts_nonme)
#[1]   376 59428
# there are 376 non-metastatic samples
dat<-rbind(counts_me,counts_nonme)
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
#[1] 59427   449
Group <- factor(rep(c("me","non_me"),times=c(73,376)),levels = c("non_me","me"))
table(Group)
#Group
#non_me     me 
#376     73
design <- model.matrix(~0+Group)
colnames(design)= levels(Group)
rownames(design)=colnames(dataexpr)
write.csv(dataexpr,"dataexpr.csv",quote = F)
#dataexpr<-read.csv("dataexpr.csv",header=T,row.names=1)

library("limma")
library("tidyverse")
library("stringr")
library("edgeR")
table(is.na(dataexpr))
#FALSE 
#26682723
dge <- DGEList(counts=dataexpr)
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
logFC_t=0
P.Value_t = 0.05
k1 = (DEG$P.Value < P.Value_t)&(DEG$logFC < -logFC_t)
k2 = (DEG$P.Value < P.Value_t)&(DEG$logFC > logFC_t)
change = ifelse(k1,"DOWN",ifelse(k2,"UP","stable"))
DEG$change <- change
DEG_limma <- DEG
#table(DEG_limma$change)
#DOWN stable     UP 
#2497  55070   1860
table(DEG_limma$change)
#DOWN stable     UP 
#6084  47725   5618 
write.csv(DEG_limma,"deg_Me_nonMe.csv",quote=F)

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
pdf("DEG_Me_nonMe_VP.pdf")
p
dev.off()

library("DOSE")
library("org.Hs.eg.db")
library("topGO")
library("clusterProfiler")
library("pathview")
###############################################up DEG functional enrichment analysis#################################################
deg_up<-subset(deg, change == "UP")
up_genes <- rownames(deg_up)
test = bitr(up_genes, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
head(test,2)
go <- enrichGO(test$ENTREZID, 
               OrgDb = org.Hs.eg.db, 
               ont='ALL',
               pAdjustMethod = 'BH',
               pvalueCutoff = 0.05, 
               qvalueCutoff = 0.05,
               keyType = 'ENTREZID')
dim(go)
write.csv(summary(go),"GO_enrich_up.csv",row.names =FALSE)
go_BP <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.05,
                  keyType = 'ENTREZID')
go_MF <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='MF',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.05,
                  keyType = 'ENTREZID')
go_CC <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='CC',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.05,
                  keyType = 'ENTREZID')
#随后对富集结果进行总览，查看BP，CC，MF的个数
dim(go[go$ONTOLOGY=='BP',])
dim(go[go$ONTOLOGY=='CC',])
dim(go[go$ONTOLOGY=='MF',])
#结果可视化
pdf("GO_BP_enrich_Me_nonMe_up.pdf")
barplot(go_BP,showCategory=10,drop=T,font.size= 12, title="GO enrichment analysis (BP)")
dev.off()
pdf("GO_MF_enrich_Me_nonMe_up.pdf")
barplot(go_MF,showCategory=10,drop=T,font.size= 12, title="GO enrichment analysis (MF)")
dev.off()
pdf("GO_CC_enrich_Me_nonMe_up.pdf")
barplot(go_CC,showCategory=10,drop=T,font.size= 12, title="GO enrichment analysis (CC)")
dev.off()
kk <- enrichKEGG(gene = test$ENTREZID,
                 organism = 'hsa', #KEGG可以用organism = 'hsa'
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
write.csv(summary(kk),"KEGG_enrich_up.csv",row.names =FALSE)
pdf("KEGG_enrich_up.pdf")
barplot(kk,showCategory=10,title="KEGG enrichment analysis")
dev.off()
#############################################################################################################################################
###############################################down DEG functional enrichment analysis#################################################
deg_down<-subset(deg, change == "DOWN")
down_genes <- rownames(deg_down)
test = bitr(down_genes, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
head(test,2)
go <- enrichGO(test$ENTREZID, 
               OrgDb = org.Hs.eg.db, 
               ont='ALL',
               pAdjustMethod = 'BH',
               pvalueCutoff = 0.05, 
               qvalueCutoff = 0.05,
               keyType = 'ENTREZID')
dim(go)
write.csv(summary(go),"GO_enrich_down.csv",row.names =FALSE)
go_BP <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.05,
                  keyType = 'ENTREZID')
go_MF <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='MF',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.05,
                  keyType = 'ENTREZID')
go_CC <- enrichGO(test$ENTREZID, 
                  OrgDb = org.Hs.eg.db, 
                  ont='CC',
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.05,
                  keyType = 'ENTREZID')
#随后对富集结果进行总览，查看BP，CC，MF的个数
dim(go[go$ONTOLOGY=='BP',])
dim(go[go$ONTOLOGY=='CC',])
dim(go[go$ONTOLOGY=='MF',])
#结果可视化
pdf("GO_BP_enrich_Me_nonMe_down.pdf")
barplot(go_BP,showCategory=10,drop=T,font.size= 12, title="GO enrichment analysis (BP)")
dev.off()
pdf("GO_MF_enrich_Me_nonMe_down.pdf")
barplot(go_MF,showCategory=10,drop=T,font.size= 12, title="GO enrichment analysis (MF)")
dev.off()
pdf("GO_CC_enrich_Me_nonMe_down.pdf")
barplot(go_CC,showCategory=10,drop=T,font.size= 12, title="GO enrichment analysis (CC)")
dev.off()
kk <- enrichKEGG(gene = test$ENTREZID,
                 organism = 'hsa', #KEGG可以用organism = 'hsa'
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
write.csv(summary(kk),"KEGG_enrich_down.csv",row.names =FALSE)
pdf("KEGG_enrich_down.pdf")
barplot(kk,showCategory=10,title="KEGG enrichment analysis")
dev.off()
#############################################################################################################################################


#The TNM Staging System
##Primary tumor (T)
#TX: Main tumor cannot be measured.
#T0: Main tumor cannot be found.
#T1, T2, T3, T4: Refers to the size and/or extent of the main tumor. The higher the number after the T, the larger the tumor or the more it has grown into nearby tissues. T's may be further divided to provide more detail, such as T3a and T3b.
##Regional lymph nodes (N)
#NX: Cancer in nearby lymph nodes cannot be measured.
#N0: There is no cancer in nearby lymph nodes.
#N1, N2, N3: Refers to the number and location of lymph nodes that contain cancer. The higher the number after the N, the more lymph nodes that contain cancer.
##Distant metastasis (M)
#MX: Metastasis cannot be measured.
#M0: Cancer has not spread to other parts of the body.
#M1: Cancer has spread to other parts of the body.



