setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/05_ImmuneLandscape/01_ssGSEA")
rm(list = ls())

# ref: https://www.jianshu.com/p/bb45958d2e9b?u_atoken=aa48a13b-57e0-4e54-99c0-12511f90c839&u_asession=015hcFMPkvQWmEqSdtYeCHiTo0c9_VfkzV_fCHgHC752PaB03g-YrY5OmclOPs9o8oX0KNBwm7Lovlpxjd_P_q4JsKWYrT3W_NKPr8w6oU7K9tLhj5qG0CcK0Hsl3CjKnhp6-k6uIK8qd-Te2lZdEIwmBkFo3NEHBv0PZUm6pbxQU&u_asig=05GMhfeFGOEHiWf9nfeNuARRkogZ0OpktbluMtPVd5FNLAbJtHH9Igs8E3vWP3hLSR0lRkCuaCT08XkrMQWesuPtAB7eXNpUE8jDxjDnfMLV5MsIYTGuP2VS0yH957Gj-tixWqTDeaCUSS_wbgyqA0qOs0yigFfZYdOgD3k98krtz9JS7q8ZD7Xtz2Ly-b0kmuyAKRFSVJkkdwVUnyHAIJzZ_xhnOgNGVDkqdKvj-YFdrHRX16bxi55mT-tAVsNXIJ9ku645VqGlBHykiFwmnEve3h9VXwMyh6PgyDIVSG1W_7V6wbSHZn9kxyd7iqk7yAH-LLY_Nj0_-S2jVka7GV0mayIZB78oqg_QAoS-wMYBqGchPoisyQRvU6rvjSOQaCmWspDxyAEEo4kbsryBKb9Q&u_aref=EscckrpHEAXPk%2Fz01TNhz7ypPyI%3D

library("GSVA")
library("limma")
library("pheatmap")
library("GSEABase")

# 读取从 GSEA 官网下载的通路数据
c2gmt <- getGmt("h.all.v2022.1.Hs.symbols.gmt")
dataexpr<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/Overlap_genes/FPKM_mRNAExp.Me_NoMe_selected.csv",
                   header=TRUE, row.names = 1)
res.ssgsea <- gsva(as.matrix(dataexpr), c2gmt, method = "ssgsea", kcdf = "Poisson", min.sz = 10)
write.csv(res.ssgsea, "ssGSEA_results.csv", quote = F)
ssgsea_t<-t(res.ssgsea)
ssgsea_t<-as.data.frame(ssgsea_t)
ssgsea_t$sample_name<-row.names(ssgsea_t)
ssgsea_t$sample_name<-substr(ssgsea_t$sample_name,1,12)
ssgsea_t$sample_name<-gsub("\\.","-",ssgsea_t$sample_name)

risk<-read.table("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/04_multivariate_Cox_regression_analysis/risk.txt",
                 sep="\t", header=T, row.names = 1)
head(risk)
risk$sample_name<-row.names(risk)
risk$sample_name<-substr(risk$sample_name,1,12)
head(risk)

merged_data<-merge(ssgsea_t, risk, by = "sample_name")
merged_data<-merged_data[,-(52:57)]
merged_data<-merged_data[,-53]
for(i in 1:nrow(merged_data)){
  merged_data$sample_name[i]<-paste(merged_data$sample_name[i], i, sep = "_")
}
rownames(merged_data)<-merged_data$sample_name
dim(merged_data)
merged_data[1:5,1:5]
merged_data<-merged_data[,-1]



library("devtools")
#devtools::install_local("D:/MyProjects/CTC/CTC_Prognosis_COAD/05_ImmuneLandscape/02_ssGSEA/Rpackage_ggcor/ggcor-master/ggcor-master",
#                        force = TRUE)
library(ggcor)
packageVersion('ggcor')
#[1] ‘0.9.4.3’
library(dplyr)
library(vegan)
library(ggplot2)
library(patchwork)

pdf("ssGSEA_hallmark_riskScore.pdf")
mantel02 <- fortify_mantel(merged_data[,51], merged_data[,1:50]) %>% 
  mutate(R = cut(r, breaks = c(-Inf, 0.08923, 0.11721, 0.13246, Inf), 
                 labels = c("<0.08923", "0.08923-0.11721", "0.11721-0.13246",">=0.13246"),
                 right = FALSE),
         p.value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">=0.05"),
                       right = FALSE))

quickcor(merged_data, type = "upper") + 
  geom_colour() + 
  scale_fill_gradientn(colours = c("blue","white" ,"red")) +
  add_link(mantel02, mapping = aes(colour = p.value),
           diag.label = TRUE)  + remove_axis("x")
dev.off()

