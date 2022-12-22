setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/07_ChemotherapyResistanceAnalysis/GSE143985")
rm(list=ls())


com_auc<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/03_Gaussian_mixture_model-based_hierarchical_clustering_method/cluster_results.csv",
                  header=T)
com_auc_sorted=com_auc[order(com_auc$auc,decreasing = T),]
genes<-com_auc_sorted$gene[1]
gene_list<-strsplit(genes,split="\\ \\+\\ ")[[1]]
uniSigExp<-read.table("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/01_univariate_Cox_regression_analysis/uniSigExp.txt",
                      sep="\t",header=T)
uniSigExp_filtered<-uniSigExp[,c("id","futime","fustat",gene_list)]
dim(uniSigExp_filtered)
#[1] 521   7
head(uniSigExp_filtered)

library("survival")
rt=uniSigExp_filtered
rownames(rt)<-rt$id
rt<-rt[,-1]
#COX模型构建
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

#输出模型参数
outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

clinical_infor<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/05_C-index/GSE143985_GPL570/clinical_infor.csv")
exprSet_filter<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/05_C-index/GSE143985_GPL570/expData_validation_1.csv",
                         row.names = 1)
genes<-c("TRIP10","NGFR","SLC48A1","SRMS")
exprSet_selected<-exprSet_filter[genes,]
exprSet_selected_t<-t(exprSet_selected)
exprSet_selected_t<-as.data.frame(exprSet_selected_t)
exprSet_selected_t$sample_name<-row.names(exprSet_selected_t)

merged_data<-merge(clinical_infor, exprSet_selected_t, by = "sample_name")
rt<-subset(merged_data, select = c("sample_name", "dfs_time", "dfs_event", "TRIP10", "NGFR", "SLC48A1", "SRMS"))
colnames(rt)<-c("id", "futime", "fustat", "TRIP10", "NGFR", "SLC48A1", "SRMS")
rownames(rt)<-rt$id
#输出病人风险值
riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
            file="risk.txt",
            sep="\t",
            quote=F,
            row.names=F)


