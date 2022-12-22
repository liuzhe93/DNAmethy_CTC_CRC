setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/04_FunctionalEnrichmentAnalysis/05_GSEAEnrichment")
rm(list = ls())

DEGs<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/04_FunctionalEnrichmentAnalysis/04_MetascapeEnrichment/deg_High_LowRisk.csv")

library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
gene<-str_trim(DEGs$X,"both") #定义gene
#开始ID转换
gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") #会有部分基因数据丢失，或者ENSEMBL
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene_df <- data.frame(logFC=DEGs$logFC, #可以是foldchange
                      SYMBOL = DEGs$X) #记住你的基因表头名字
gene_df <- merge(gene_df,gene,by="SYMBOL")
head(gene_df)

geneList<-gene_df$logFC #第二列可以是folodchange，也可以是logFC
names(geneList)=gene_df$ENTREZID #使用转换好的ID
geneList=sort(geneList,decreasing = T) #从高到低排序

kegmt<-read.gmt("c2.cp.kegg.v2022.1.Hs.entrez.gmt") #读gmt文件
KEGG<-GSEA(geneList,TERM2GENE = kegmt,pvalueCutoff = 1) #GSEA分析
dim(KEGG)
write.csv(KEGG, "KEGG_results.csv", quote = F)
library(ggplot2)
#dotplot(KEGG) #出点图 
#dotplot(KEGG,color="pvalue")  #按p值出点图 
#dotplot(KEGG,split=".sign")+facet_grid(~.sign) #出点图，并且分面激活和抑制
#dotplot(KEGG,split=".sign")+facet_wrap(~.sign,scales = "free") #换个显示方式
library(enrichplot)
#特定通路作图
#gseaplot2(KEGG,1,color="red",pvalue_table = T) # 按第一个做二维码图，并显示p值
#gseaplot2(KEGG,1:10,color="red") #按第一到第十个出图，不显示p值
gseaplot2(KEGG,2,color="red",pvalue_table = T)


pdf("GSEA_highriskGroup.pdf")
gseaplot2(KEGG,2:6,color="red",pvalue_table = T)
dev.off()

pdf("GSEA_lowriskGroup.pdf")
gseaplot2(KEGG,c(1,23,25,26,35),color="red",pvalue_table = T)
dev.off()


