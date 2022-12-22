setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/01_univariate_Cox_regression_analysis")
rm(list=ls())

pFilter=0.05                                                      #定义单因素显著性
library(survival)                                                 #引用包
outTab<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/Overlap_genes/corResult.csv",header=T)        #输出相关性结果
genes<-outTab$FPKM
FPKM_TCGA<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/Data_Collection/TCGA-COAD/mRNA_exp/FPKM.csv")
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
FPKM_anno<-merge(anno,FPKM_TCGA,by.x="gene_id",by.y="X")
FPKM_anno<-FPKM_anno[,-1]
FPKM_anno<-aggregate(FPKM_anno,by=list(FPKM_anno$gene_name),FUN=median)
write.csv(FPKM_anno,"FPKM_anno_uniq.csv",quote=F)
#FPKM_anno<-read.csv("FPKM_anno_uniq.csv",header=T,row.names=1)
dim(FPKM_anno)
#[1] 59427   521
FPKM_anno_t<-t(FPKM_anno)
FPKM_anno_t<-as.data.frame(FPKM_anno_t)
FPKM_selected<-FPKM_anno_t[,genes]
FPKM_selected$sample_name<-gsub("\\.","-",rownames(FPKM_selected))
FPKM_selected$sample_name<-substr(FPKM_selected$sample_name,1,12)
result_time<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/Metastatic_nonMetastatic/result_time.csv",header=T)
exptime<-merge(FPKM_selected,result_time,by.x="sample_name",by.y="submitter_id")
dim(exptime)
#[1] 521  173
head(exptime)
expTime_data<-cbind(exptime[,c("sample_name","futime","fustat")],exptime[,2:167])
for(i in 1:521){
  rownames(expTime_data)[i]<-paste(expTime_data$sample_name[i],i,sep="_")
}
expTime_data<-expTime_data[,-1]
head(expTime_data)
outTab=data.frame()
rt<-expTime_data
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     z=coxSummary$coefficients[,"z"],
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
outTab = outTab[is.na(outTab$pvalue)==FALSE,]
outTab=outTab[order(as.numeric(as.vector(outTab$pvalue))),]
write.table(outTab,file="uniCoxResult.txt",sep="\t",row.names=F,quote=F)
sigTab=outTab[as.numeric(as.vector(outTab$pvalue))<pFilter,]
write.table(sigTab,file="uniCoxResult.Sig.txt",sep="\t",row.names=F,quote=F)
dim(outTab)
#[1] 166   6
dim(sigTab)
#[1] 25  6
sigGenes=c("futime","fustat")
sigGenes=c(sigGenes,as.vector(sigTab[,1]))
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)
uniSigExp<-read.table("uniSigExp.txt",sep="\t",header=T)
dim(uniSigExp)
#[1] 521  28

library("ggplot2")
head(outTab)
outTab$z<-as.numeric(outTab$z)
outTab$pvalue<-as.numeric(outTab$pvalue)
outTab$color<-ifelse(outTab$pvalue<0.05 ,"red","grey")
color<-c(red = "red", grey = "grey")
p<-ggplot(outTab, aes(z, -log10(pvalue), col = color)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = color) +
  labs(x="Univariate Cox coefficient", y = "-log10(P-Value)") +
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
pdf("UnivariateCOXRegressionAnalysis.pdf")
p
dev.off()
