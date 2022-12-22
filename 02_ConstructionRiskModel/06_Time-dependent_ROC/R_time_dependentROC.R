setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/06_Time-dependent_ROC")
rm(list=ls())

library(timeROC)
library(survival)
library(survivalROC)
###########################################TCGA数据的ROC曲线分析####################################################
risk<-read.table("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/04_multivariate_Cox_regression_analysis/risk.txt",
                      sep="\t",header=T)
dim(risk)
#[1] 521   9
head(risk)
time_roc_res <- timeROC(
  T = risk$futime,
  delta = risk$fustat,
  marker = risk$riskScore,
  cause = 1,
  weighting="marginal",
  times = c(3, 5, 10),
#  times = c(1,2,3,4,5,6,7,8,9,10),
  ROC = TRUE,
  iid = TRUE
)
#time_roc_res$AUC
#     t=1       t=2       t=3       t=4       t=5       t=6       t=7       t=8       t=9      t=10 
#0.7279297 0.7024452 0.6824516 0.6822810 0.6389583 0.5999614 0.6032997 0.5819867 0.6345913 0.6221901
time_roc_res$AUC
#     t=3       t=5      t=10 
#0.6824516 0.6389583 0.6221901 
confint(time_roc_res, level = 0.95)$CI_AUC
#      2.5% 97.5%
#t=3  60.52 75.97
#t=5  53.56 74.23
#t=10 41.99 82.45

time_ROC_df <- data.frame(
  TP_3year = time_roc_res$TP[, 1],
  FP_3year = time_roc_res$FP[, 1],
  TP_5year = time_roc_res$TP[, 2],
  FP_5year = time_roc_res$FP[, 2],
  TP_10year = time_roc_res$TP[, 3],
  FP_10year = time_roc_res$FP[, 3]
)
library(ggplot2)
pdf("ROC_curve_year3510_TCGA.pdf")
ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_10year, y = TP_10year), size = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 10 years = ", sprintf("%.3f", time_roc_res$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(face = "bold", size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
  )
dev.off()

###############################validation dataset 1的ROC曲线分析####################################################
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

rt=mydata_vali_1
rownames(rt)<-rt$sample_name
rt<-rt[,-1]
#COX模型构建
multiCox=coxph(Surv(dfs_time, dfs_event) ~ ., data = rt)
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
write.table(outTab,file="multiCox_GSE143985.xls",sep="\t",row.names=F,quote=F)
riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("dfs_time","dfs_event",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
            file="risk_GSE143985.txt",
            sep="\t",
            quote=F,
            row.names=F)
exp_dataset1<-read.table("risk_GSE143985.txt",header = T, row.names = 1)
time_roc_res2 <- timeROC(
  T = exp_dataset1$dfs_time,
  delta = exp_dataset1$dfs_event,
  marker = exp_dataset1$riskScore,
  cause = 1,
  weighting="marginal",
  times = c(3, 5, 10),
#  times = c(1,2,3,4,5,6,7,8,9,10),
  ROC = TRUE,
  iid = TRUE
)
#time_roc_res2$AUC
#     t=1       t=2       t=3       t=4       t=5       t=6       t=7       t=8       t=9      t=10 
#0.8105810 0.7692289 0.8180662 0.8127152 0.8124294 0.7965145 0.7917701 0.7548633 0.7313885 0.7069693
time_roc_res2$AUC
#     t=3       t=5      t=10 
#0.8180662 0.8124294 0.706969
confint(time_roc_res2, level = 0.95)$CI_AUC
#      2.5% 97.5%
#t=3  71.53 92.09
#t=5  70.78 91.71
#t=10 45.86 95.53
time_ROC_df <- data.frame(
  TP_3year = time_roc_res2$TP[, 1],
  FP_3year = time_roc_res2$FP[, 1],
  TP_5year = time_roc_res2$TP[, 2],
  FP_5year = time_roc_res2$FP[, 2],
  TP_10year = time_roc_res2$TP[, 3],
  FP_10year = time_roc_res2$FP[, 3]
)
library(ggplot2)
pdf("ROC_curve_year3510_GSE143985.pdf")
ggplot(data = time_ROC_df) +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_10year, y = TP_10year), size = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
           x = 0.75, y = 0.25, size = 4.5,
           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res2$AUC[[1]])), color = "#BC3C29FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.15, size = 4.5,
           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res2$AUC[[2]])), color = "#0072B5FF"
  ) +
  annotate("text",
           x = 0.75, y = 0.05, size = 4.5,
           label = paste0("AUC at 10 years = ", sprintf("%.3f", time_roc_res2$AUC[[3]])), color = "#E18727FF"
  ) +
  labs(x = "False positive rate", y = "True positive rate") +
  theme(
    axis.text = element_text(face = "bold", size = 11, color = "black"),
    axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
    axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
  )
dev.off()
###############################validation dataset 2的ROC曲线分析####################################################
#clinical_infor<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/05_C-index/GSE63624_GPL5175/clinical_infor.csv",
#                         header = T)
#data.exprs<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/05_C-index/GSE63624_GPL5175/expData_validation_2.csv",
#                         header = T, row.names = 1)
#com_auc<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/03_Gaussian_mixture_model-based_hierarchical_clustering_method/cluster_results.csv",
#                  header=T)
#com_auc_sorted=com_auc[order(com_auc$auc,decreasing = T),]
#genes<-com_auc_sorted$gene[1]
#gene_list<-strsplit(genes,split="\\ \\+\\ ")[[1]]

#TRIP10_exp<-data.exprs["3818515",]
#NGFR_exp<-data.exprs["3725685",]
#SLC48A1_exp<-data.exprs["3413212",]
#SRMS_exp<-data.exprs["3913945",]
#exp_selected<-rbind(TRIP10_exp,NGFR_exp,SLC48A1_exp,SRMS_exp)
#exp_selected_t<-t(exp_selected)
#exp_selected_t<-as.data.frame(exp_selected_t)
#colnames(exp_selected_t)<-c("TRIP10","NGFR","SLC48A1","SRMS")
#exp_selected_t$sample_name<-row.names(exp_selected_t)
#mydata_vali_2<-merge(clinical_infor,exp_selected_t,by="sample_name")
#mydata_vali_2$rfs_event<-as.integer(mydata_vali_2$rfs_event)
#row.names(mydata_vali_2)<-mydata_vali_2$sample_name
#mydata_vali_2<-mydata_vali_2[,-1]
#rt=mydata_vali_2
#COX模型构建
#multiCox=coxph(Surv(rfs_time, rfs_event) ~ ., data = rt)
#multiCox=step(multiCox,direction = "both")
#multiCoxSum=summary(multiCox)
#输出模型参数
#outTab=data.frame()
#outTab=cbind(
#  coef=multiCoxSum$coefficients[,"coef"],
#  HR=multiCoxSum$conf.int[,"exp(coef)"],
#  HR.95L=multiCoxSum$conf.int[,"lower .95"],
#  HR.95H=multiCoxSum$conf.int[,"upper .95"],
#  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
#outTab=cbind(id=row.names(outTab),outTab)
#outTab=gsub("`","",outTab)
#write.table(outTab,file="multiCox_GSE63624.xls",sep="\t",row.names=F,quote=F)
#riskScore=predict(multiCox,type="risk",newdata=rt)
#coxGene=rownames(multiCoxSum$coefficients)
#coxGene=gsub("`","",coxGene)
#outCol=c("rfs_time","rfs_event",coxGene)
#risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
#write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
#            file="risk_GSE63624.txt",
#            sep="\t",
#            quote=F,
#            row.names=F)
#exp_dataset2<-read.table("risk_GSE63624.txt",header = T, row.names = 1)
#time_roc_res3 <- timeROC(
#  T = exp_dataset2$rfs_time,
#  delta = exp_dataset2$rfs_event,
#  marker = exp_dataset2$riskScore,
#  cause = 1,
#  weighting="marginal",
#  times = c(1, 3, 5),
#  times = c(1,2,3,4,5,6),
#  ROC = TRUE,
#  iid = TRUE
#)
#time_roc_res3$AUC
#      t=1       t=2       t=3       t=4       t=5       t=6 
#0.5322581 0.6405219 0.6777238 0.6550374 0.6384537 0.7083489
#time_roc_res3$AUC
#      t=1       t=3       t=5 
#0.5322581 0.6777238 0.6384537
confint(time_roc_res3, level = 0.95)$CI_AUC
#     2.5% 97.5%
#t=1 33.85 72.60
#t=3 51.83 83.72
#t=5 45.38 82.31
#time_ROC_df <- data.frame(
#  TP_1year = time_roc_res3$TP[, 1],
#  FP_1year = time_roc_res3$FP[, 1],
#  TP_3year = time_roc_res3$TP[, 2],
#  FP_3year = time_roc_res3$FP[, 2],
#  TP_5year = time_roc_res3$TP[, 3],
#  FP_5year = time_roc_res3$FP[, 3]
#)
#library(ggplot2)
#pdf("ROC_curve_year135_GSE63624.pdf")
#ggplot(data = time_ROC_df) +
#  geom_line(aes(x = FP_1year, y = TP_1year), size = 1, color = "#BC3C29FF") +
#  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#0072B5FF") +
#  geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#E18727FF") +
#  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
#  theme_bw() +
#  annotate("text",
#           x = 0.75, y = 0.25, size = 4.5,
#           label = paste0("AUC at 1 years = ", sprintf("%.3f", time_roc_res3$AUC[[1]])), color = "#BC3C29FF"
#  ) +
#  annotate("text",
#           x = 0.75, y = 0.15, size = 4.5,
#           label = paste0("AUC at 3 years = ", sprintf("%.3f", time_roc_res3$AUC[[2]])), color = "#0072B5FF"
#  ) +
#  annotate("text",
#           x = 0.75, y = 0.05, size = 4.5,
#           label = paste0("AUC at 5 years = ", sprintf("%.3f", time_roc_res3$AUC[[3]])), color = "#E18727FF"
#  ) +
#  labs(x = "False positive rate", y = "True positive rate") +
#  theme(
#    axis.text = element_text(face = "bold", size = 11, color = "black"),
#    axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
#    axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
#  )
#dev.off()
###############################比较两组time-dependent AUC####################################################
mydata<-cbind(time_roc_res$AUC,time_roc_res2$AUC)
colnames(mydata)<-c("TCGA","GSE143985")
mydata<-as.data.frame(mydata)
mydata$Year<-c(1,2,3,4,5,6,7,8,9,10)
all_data<-reshape2::melt(mydata, id.vars = "Year", variable.names = "DataSet")
write.csv(all_data, "AUC_timedependent.csv", row.names = F, quote = F)
#all_data<-read.csv("AUC_timedependent.csv",header=T)
pdf("tROC_2DatasetMerged.pdf")
theme_set(theme_bw())
p<-ggplot(data = all_data, mapping = aes(x = Year, y = value, colour = variable)) + geom_line() + ylim(c(0,1))
p_bottom=p+theme(legend.position = "top")
p_bottom
dev.off()







