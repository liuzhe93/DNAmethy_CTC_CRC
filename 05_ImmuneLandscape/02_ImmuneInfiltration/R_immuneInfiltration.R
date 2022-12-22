setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/05_ImmuneLandscape/02_ImmuneInfiltration")
rm(list = ls())
#参考资料：https://www.jianshu.com/p/03a7440c0960

library(tinyarray)
library(tidyverse)

#######################临床数据的处理###################################
result_time<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/Metastatic_nonMetastatic/result_time.csv",header=T)
head(result_time)
result_time$Age_ope<-ifelse(result_time$age_at_index>=60,">=60","<60")
result_time$Stage_ope<-ifelse(result_time$ajcc_pathologic_stage=="Stage IA", "Stage I",
                              ifelse(result_time$ajcc_pathologic_stage=="Stage IIA", "Stage II",
                                     ifelse(result_time$ajcc_pathologic_stage=="Stage IIB", "Stage II",
                                            ifelse(result_time$ajcc_pathologic_stage=="Stage IIC", "Stage II",
                                                   ifelse(result_time$ajcc_pathologic_stage=="Stage IIIA", "Stage III",
                                                          ifelse(result_time$ajcc_pathologic_stage=="Stage IIIB", "Stage III",
                                                                 ifelse(result_time$ajcc_pathologic_stage=="Stage IIIC", "Stage III",
                                                                        ifelse(result_time$ajcc_pathologic_stage=="Stage IVA", "Stage IV",
                                                                               ifelse(result_time$ajcc_pathologic_stage=="Stage IVB", "Stage IV",
                                                                                      ifelse(result_time$ajcc_pathologic_stage=="Stage I", "Stage I",
                                                                                             ifelse(result_time$ajcc_pathologic_stage=="Stage II", "Stage II",
                                                                                                    ifelse(result_time$ajcc_pathologic_stage=="Stage III", "Stage III",
                                                                                                           ifelse(result_time$ajcc_pathologic_stage=="Stage IV", "Stage IV",
                                                                                                                  "NA")))))))))))))
result_time$Status_ope<-ifelse(result_time$fustat==0,"Alive","Dead")
result_time<-subset(result_time, select = c("submitter_id","Age_ope","Stage_ope","Status_ope"))
colnames(result_time)<-c("sample_name", "Age", "Stage", "Status")
######################risk以及risk score###################################
risk<-read.table("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/04_multivariate_Cox_regression_analysis/risk.txt",
                 sep="\t", header=T, row.names = 1)
head(risk)
risk$sample_name<-row.names(risk)
risk$sample_name<-substr(risk$sample_name,1,12)
risk<-subset(risk, select = c(riskScore, risk, sample_name))


dataexpr<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/Overlap_genes/FPKM_mRNAExp.Me_NoMe_selected.csv",
                   header=TRUE, row.names = 1)

exp2 = dataexpr
exp2 = rownames_to_column(exp2)
write.table(exp2,file = "exp.txt",row.names = F,quote = F,sep = "\t")

source("CIBERSORT.R")

if(T){
  TME.results = CIBERSORT("LM22.txt", 
                          "exp.txt" , 
                          perm = 1000, 
                          QN = F)
  save(TME.results,file = "ciber_COAD.Rdata")
}
load("ciber_COAD.Rdata")
TME.results[1:4,1:4]

re <- TME.results[,-(23:25)]
library(pheatmap)

re_df<-as.data.frame(re)
re_df$sample_name<-row.names(re_df)
re_df$sample_name<-substr(re_df$sample_name,1,12)
re_df$sample_name<-gsub("\\.","-",re_df$sample_name)
merged_re<-merge(re_df,risk,by="sample_name")
merged_re<-merge(merged_re, result_time, by = "sample_name")
merged_re_sort<-merged_re[order(merged_re[,25],decreasing=TRUE),]
merged_re_heatmap<-merged_re_sort[,-c(1,24,25,26,27,28)]
#那些在一半以上样本里丰度为0的免疫细胞，就不展示在热图里了
k <- apply(merged_re_heatmap,2,function(x) {sum(x == 0) < nrow(TME.results)/2})
table(k)
#k
#FALSE  TRUE 
#   9    13 
#re2 <- as.data.frame(t(merged_re_heatmap[,k]))
re2 <- as.data.frame(t(merged_re_heatmap))
Group<-merged_re_sort$risk
Age<-merged_re_sort$Age
Stage<-merged_re_sort$Stage
Status<-merged_re_sort$Status
RiskScore<-merged_re_sort$riskScore


an = data.frame(group = Group,
                age = Age,
                stage = Stage,
                status = Status,
                riskScore = RiskScore,
                row.names = colnames(re2))
pdf("Heatmap_immune_infiltration.pdf")
pheatmap(re2,scale = "row",
         show_colnames = F,
         annotation_col = an,cluster_rows = T,cluster_cols = F,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()


library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))

dat <- re %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)

dat$sample_name<-dat$Sample
dat$sample_name<-substr(dat$sample_name,1,12)
dat$sample_name<-gsub("\\.","-",dat$sample_name)
merged_dat<-merge(dat, risk, by = "sample_name")
merged_dat$Group = merged_dat$risk
merged_dat$Group<-ifelse(merged_dat$Group=="high","2_high","1_low")

library(ggpubr)

merged_dat_Tcellgammadelta<-subset(merged_dat,Cell_type=="T cells gamma delta")
merged_dat_Tcellgammadelta_low<-subset(merged_dat_Tcellgammadelta, risk == "low")
merged_dat_Tcellgammadelta_high<-subset(merged_dat_Tcellgammadelta, risk == "high")
mean(merged_dat_Tcellgammadelta_high$Proportion)
#[1] 0.0005259607
mean(merged_dat_Tcellgammadelta_low$Proportion)
#[1] 0.001376141

pdf("Boxplot_immune_infiltration.pdf", height = 6, width = 10)
ggplot(merged_dat,aes(Cell_type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test")
dev.off()



