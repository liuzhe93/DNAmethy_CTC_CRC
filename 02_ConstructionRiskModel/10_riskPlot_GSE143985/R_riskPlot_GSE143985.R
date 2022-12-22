setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/10_riskPlot_GSE143985/")
rm(list=ls())

clinical_infor<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/05_C-index/GSE143985_GPL570/clinical_infor.csv",header = T)
exprSet_filter<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/05_C-index/GSE143985_GPL570/expData_validation_1.csv",
                         header = T, row.names = 1)
com_auc<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/03_Gaussian_mixture_model-based_hierarchical_clustering_method/cluster_results.csv",
                  header=T)
com_auc_sorted=com_auc[order(com_auc$auc,decreasing = T),]
genes<-com_auc_sorted$gene[1]
gene_list<-strsplit(genes,split="\\ \\+\\ ")[[1]]
exp_dataset1<-exprSet_filter[gene_list,]
exp_selected_t<-t(exp_dataset1)
exp_selected_t<-as.data.frame(exp_selected_t)
exp_selected_t$sample_name<-row.names(exp_selected_t)
mydata_vali_1<-merge(clinical_infor,exp_selected_t,by="sample_name")
mydata_vali_1$dfs_event<-as.integer(mydata_vali_1$dfs_event)
risk<-read.table("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/06_Time-dependent_ROC/risk_GSE143985.txt",
                 sep="\t",header=T, row.names = 1)
risk$sample_name<-rownames(risk)
gse_data<-merge(mydata_vali_1,risk,by="sample_name")
gse_data<-subset(gse_data,select=c("sample_name","dfs_event.x","dfs_time.x","TRIP10","NGFR","SLC48A1.x","SRMS","riskScore",
                                   "risk"))
colnames(gse_data)<-c("sample_name","dfs_event","dfs_time","TRIP10","NGFR","SLC48A1","SRMS","riskScore","risk")
rownames(gse_data)<-gse_data$sample_name
gse_data<-gse_data[,-1]
  
library(pheatmap)
rt<-gse_data
rt=rt[order(rt$riskScore),]                                     #按照riskScore对样品排序

#绘制风险曲线
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
#line[line>10]=10
pdf(file="riskScore.pdf",width = 12,height = 5)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("green",lowLength),
           rep("red",highLength)))
abline(h=median(rt$riskScore),v=lowLength,lty=2)
dev.off()

#绘制生存状态图
color=as.vector(rt$dfs_event)
color[color==1]="red"
color[color==0]="green"
pdf(file="survStat.pdf",width = 12,height = 5)
plot(rt$dfs_time,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
abline(v=lowLength,lty=2)
dev.off()

#绘制风险热图
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pdf(file="heatmap.pdf",width = 12,height = 5)
pheatmap(rt1, 
         annotation=annotation, 
         cluster_cols = FALSE,
         fontsize_row=11,
         fontsize_col=3,
         color = colorRampPalette(c("blue", "white", "red"))(50) )
dev.off()

