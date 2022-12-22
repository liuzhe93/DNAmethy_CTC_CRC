setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/05_ImmuneLandscape/03_ImmuneCheckpointsExp")
rm(list = ls())

#genelist<-c("CD274","CD8A","CTLA4","CXCL10","CXCL9","GZMA","GZMB","HAVCR2","IDO1","IFNG","LAG3","PDCD1","PRF1","TBX2","TNF")
genelist<-c("CD274","CTLA4","CXCL10","GZMA","GZMB","IFNG","LAG3","PDCD1","PRF1","TBX2","TNF")

dataexpr<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/Overlap_genes/FPKM_mRNAExp.Me_NoMe_selected.csv",
                   header=TRUE, row.names = 1)
dataexpr_selected<-dataexpr[genelist,]
dataexpr_selected_t<-t(dataexpr_selected)
dataexpr_selected_t<-as.data.frame(dataexpr_selected_t)
dataexpr_selected_t$sample_name<-row.names(dataexpr_selected_t)
dataexpr_selected_t$sample_name<-substr(dataexpr_selected_t$sample_name,1,12)
dataexpr_selected_t$sample_name<-gsub("\\.","-",dataexpr_selected_t$sample_name)
head(dataexpr_selected_t)

risk<-read.table("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/04_multivariate_Cox_regression_analysis/risk.txt",
                 sep="\t", header=T, row.names = 1)
risk$sample_name<-row.names(risk)
risk$sample_name<-substr(risk$sample_name,1,12)
risk<-subset(risk, select = c(riskScore, risk, sample_name))
head(risk)

merged_data<-merge(dataexpr_selected_t, risk, by = "sample_name")
for(i in 1:nrow(merged_data)){
  merged_data$sample_name[i]<-paste(merged_data$sample_name[i],i,sep="_")
}
rownames(merged_data)<-merged_data$sample_name
merged_data<-merged_data[,-1]
merged_data<-merged_data[,-12]

dat <- merged_data %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Gene,value = FPKM,-Sample)
dat$sample_name<-substr(dat$Sample,1,12)
merged<-merge(dat,risk,by="sample_name")
merged$FPKM<-as.numeric(merged$FPKM)
merged$LOG2FPKM<-log2(merged$FPKM)
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))

merged<-subset(merged, Gene!="risk")
merged$risk<-ifelse(merged$risk=="high","2_high","1_low")



pdf("Boxplot_ImmuneCheckPoints.pdf")
ggplot(merged,aes(Gene,LOG2FPKM,fill = risk)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Gene", y = "Log 2 transformed FPKM") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ stat_compare_means(aes(group = risk,label = ..p.signif..),method = "kruskal.test")
dev.off()

