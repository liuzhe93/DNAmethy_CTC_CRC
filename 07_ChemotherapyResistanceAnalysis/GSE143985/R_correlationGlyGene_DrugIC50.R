setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/07_ChemotherapyResistanceAnalysis/GSE143985/")
rm(list=ls())

library("WGCNA")

##载入数据

ic50<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/07_ChemotherapyResistanceAnalysis/GSE143985_drugSensitivity.csv")
dataFilt_LIHC_final<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/05_C-index/GSE143985_GPL570/expData_validation_1.csv",
                              header = T, check.names = FALSE)

# 定义行名
rownames(dataFilt_LIHC_final) <- dataFilt_LIHC_final[,1]
dataFilt_LIHC_final <- dataFilt_LIHC_final[,-1]
# 先看一下矩阵长啥样，心里有个数：每一行是一个基因，每一列是一个样本
#View(dataFilt_LIHC_final)
dim(dataFilt_LIHC_final)
# [1] 20857    91

gly_genes<-c("TRIP10", "NGFR", "SLC48A1", "SRMS") 
FPKM<-dataFilt_LIHC_final[gly_genes,]
FPKM_rmNA<-na.omit(FPKM)
dim(FPKM_rmNA)
#[1]  4 91
FPKM_rmNA<-subset(FPKM,subset=(rowSums(is.na(FPKM)) == 0))
dim(FPKM_rmNA)
#[1]  4 449

merged<-t(FPKM_rmNA)
merged<-as.data.frame(merged)
row.names(merged)<-gsub("-",".",row.names(merged))
merged$SampleId<-row.names(merged)

ic50$SampleId<-ic50$X

data_merged<-merge(merged,ic50,by="SampleId")
data_input<-data_merged
row.names(data_input)<-data_merged$SampleId
data_input<-data_input[,-1]

geneExp<-data_input[,1:4]
IC50_score<-data_input[,6:65]
correlation_pvalue<-corAndPvalue(geneExp,IC50_score,method="pearson",use="p")
result_cor<-correlation_pvalue$cor
result_p<-correlation_pvalue$p

result_cor_t<-t(result_cor)
result_p_t<-t(result_p)

library(reshape2)
library(ggplot2)          
#绘制矩阵气泡图 matrix bubble plot
#注释：package使用之前需要调用
result_cor_transform<-melt(result_cor_t)
result_p_transform<-melt(result_p_t)
result_p_transform$size=-log10(result_p_transform$value)
colnames(result_cor_transform)<-c("Drug","RiskGene","correlation")
colnames(result_p_transform)<-c("Drug","RiskGene","pvalue","minuslog10_p_value")
data_processed<-cbind(result_cor_transform,result_p_transform)
data_processed<-subset(data_processed,select=c("Drug","RiskGene","correlation","minuslog10_p_value"))
pdf("bubble.pdf",height=10,width=8)
p<-ggplot(data_processed, aes(x = RiskGene, y = Drug, size = minuslog10_p_value, color=correlation)) + geom_point()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray"),
        panel.border = element_rect(colour="black",fill=NA))
p + scale_color_gradient2(low = "blue", mid = "white", high =  "red") +theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = .5))
dev.off()

-log10(0.05)
#[1] 1.30103
data_sig_byP<-data_processed[which(data_processed$minuslog10_p_value > 1.30103),]
data_sig_byCor_large<-data_sig_byP[which(data_sig_byP$correlation > 0.4),]
data_sig_byCor_less<-data_sig_byP[which(data_sig_byP$correlation < (-0.4)),]




