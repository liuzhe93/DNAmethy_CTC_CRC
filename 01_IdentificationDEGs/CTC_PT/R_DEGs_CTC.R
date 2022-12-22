setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/CTC_PT")
remove(list=ls())

data_FPKM<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/Data_Collection/CTC_primaryTumor/GSE74369_COMPLETE_RAW_FPKM_TABLE.csv",
                    row.names = 1, header = T)
dim(data_FPKM)
#[1] 21921    53
# there are 21921 genes and 53 samples
data_label<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/Data_Collection/CTC_primaryTumor/SampleLabel.csv", header = T)
table(data_label$SampleLabel)
#PrimaryTumor     TECC/CTM 
#9           18
# there are 9 primary tumor samples and 18 TECC/CTM samples
data_FPKM_selected_samples<-subset(data_FPKM, select = data_label$SampleName)
data_FPKM_selected_samples<-data.frame(data_FPKM_selected_samples[,19:27],data_FPKM_selected_samples[,1:18])
write.csv(data_FPKM_selected_samples,"FPKM_CTC_PT.csv",quote=F)

library("limma")
expMatrix <- data_FPKM_selected_samples
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpms <- apply(expMatrix,2,fpkmToTpm)
tpms[1:3,]
colSums(tpms)
group_list=c(rep('PT',9),rep('CTC',18))
group_list <- factor(group_list,levels = c("PT","CTC"),ordered = F)
exprSet <- tpms
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
exprSet <- log2(exprSet+1)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
dat <- exprSet
design=model.matrix(~factor( group_list ))
fit=lmFit(dat,design)
fit=eBayes(fit)
options(digits = 4)
topTable(fit,coef=2,adjust='BH')
bp=function(g){
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
deg=topTable(fit,coef=2,adjust='BH',number = Inf)
head(deg) 
save(deg,file = 'deg.Rdata')
#logFC_t=0.5849625
#deg$g=ifelse(deg$P.Value>0.05,'stable',
#             ifelse( deg$logFC > logFC_t,'UP',
#                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
#)
#table(deg$g)
#head(deg)
#table(deg$g)
#DOWN stable     UP 
#2744  17592   1585
logFC_t=0
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)
head(deg)
table(deg$g)
#DOWN stable     UP 
#2822  12558   6541
write.csv(deg,"deg_CTC_PT.csv",quote = F)

library("ggplot2")
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
pdf("DEG_CTC_PT_VP.pdf")
p
dev.off()

library("DOSE")
library("org.Hs.eg.db")
library("topGO")
library("clusterProfiler")
library("pathview")

###############################################up DEG functional enrichment analysis#################################################
deg_up<-subset(deg, g == "UP")
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
#[1]  0 10
dim(go[go$ONTOLOGY=='CC',])
#[1]  0 10
dim(go[go$ONTOLOGY=='MF',])
#[1]  6 10
#结果可视化
pdf("GO_BP_enrich_CTC_PT_up.pdf")
barplot(go_BP,showCategory=10,drop=T,font.size= 12, title="GO enrichment analysis (BP)")
dev.off()
pdf("GO_MF_enrich_CTC_PT_up.pdf")
barplot(go_MF,showCategory=10,drop=T,font.size= 12, title="GO enrichment analysis (MF)")
dev.off()
pdf("GO_CC_enrich_CTC_PT_up.pdf")
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
deg_down<-subset(deg, g == "DOWN")
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
#[1] 787  10
dim(go[go$ONTOLOGY=='CC',])
#[1] 204  10
dim(go[go$ONTOLOGY=='MF',])
#[1] 78 10
#结果可视化
pdf("GO_BP_enrich_CTC_PT_down.pdf")
barplot(go_BP,showCategory=10,drop=T,font.size= 12, title="GO enrichment analysis (BP)")
dev.off()

pdf("GO_MF_enrich_CTC_PT_down.pdf")
barplot(go_MF,showCategory=10,drop=T,font.size= 12, title="GO enrichment analysis (MF)")
dev.off()

pdf("GO_CC_enrich_CTC_PT_down.pdf")
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



