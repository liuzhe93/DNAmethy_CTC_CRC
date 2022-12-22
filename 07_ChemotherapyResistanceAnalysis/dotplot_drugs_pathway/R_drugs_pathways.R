setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/07_ChemotherapyResistanceAnalysis/dotplot_drugs_pathway")
rm(list = ls())

query_result<-read.table("D:/MyProjects/CTC/CTC_Prognosis_COAD/07_ChemotherapyResistanceAnalysis/CMap/DEGs_High_LowGroup/CMap/TCGA_COAD_202210122159/my_analysis.sig_queryl1k_tool.6346c843ad37c70012bc03d5/my_analysis.sig_queryl1k_tool.6346c843ad37c70012bc03d5/arfs/TAG/query_result.gct",
                         skip = 2, sep = "\t", header = T)
query_result<-query_result[-1,]
dim(query_result)
#[1] 59308    21
query_result<-subset(query_result, pert_type == "trt_cp")
dim(query_result)
#[1] 20631    21
query_result<-subset(query_result, moa != "-666")
dim(query_result)
#[1] 7515   21
query_result<-subset(query_result, target_name != "-666")
dim(query_result)
#[1] 5866   21
#################这里，我们根据norm_cs排序，选择排名前60的用于后续分析#####################################
head60<-head(query_result,60)
length(unique(head60$moa))
#[1] 50
length(unique(head60$pert_iname))
#[1] 58

mydata<-subset(head60, select = c("pert_iname", "moa"))
write.csv(unique(mydata$moa),"moa.csv")
write.csv(unique(mydata$pert_iname),"pert_iname.csv")
t1<-unique(mydata$pert_iname) #构建行名
t2<-unique(mydata$moa) #构建列名
length(t1)
#[1] 58
length(t2)
#[1] 50
#58*50
datamatrix<-matrix(0,58,50)
rownames(datamatrix)<-t1
colnames(datamatrix)<-t2
for (i in 1:nrow(mydata)) {
  datamatrix[mydata[i,1],mydata[i,2]]<-1 #对应位置打上1
}
write.table(datamatrix,"datamatrix.txt",sep = "\t",quote = FALSE)
barplot_data<-colSums(datamatrix)
barplot_data<-as.data.frame(barplot_data)
barplot_data$moa<-row.names(barplot_data)
barplot_data$counts<-barplot_data$barplot_data
barplot_data<-subset(barplot_data, select = c("moa", "counts"))
barplot_data<-barplot_data[order(barplot_data$counts, decreasing = T),]
write.csv(barplot_data, "Barplot_dat.csv", row.names = F, quote = F)
barplot_data<-barplot_data[order(barplot_data$counts, decreasing = F),]
pdf("barplot_MOAcounts.pdf")
barplot(barplot_data$counts)
dev.off()

merged_data<-merge(mydata, barplot_data, by = "moa")
merged_data_sort<-merged_data[order(merged_data$counts, decreasing = T),]

merged_data_sort$moa<-factor(merged_data_sort$moa, levels = barplot_data$moa)

library(ggplot2)
pdf("bubble_sorted.pdf", height = 10, width = 20)
p<-ggplot(merged_data_sort, aes(x = pert_iname, y = moa, size = 0.01)) + geom_point()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray"),
        panel.border = element_rect(colour="black",fill=NA))
p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = .5))
dev.off()















###########################################
query_result_HDAC_inhibitor<-subset(query_result, moa =="HDAC inhibitor")
query_result_Topoisomerase_inhibitor<-subset(query_result, moa =="Topoisomerase inhibitor")
query_result_CDK_inhibitor<-subset(query_result, moa =="CDK inhibitor")
query_result_JNK_inhibitor<-subset(query_result, moa =="JNK inhibitor")
query_result_MTOR_inhibitor<-subset(query_result, moa =="MTOR inhibitor")
query_result_RNA_polymerase_inhibitor<-subset(query_result, moa =="RNA polymerase inhibitor")
query_result_Tyrosine_kinase_inhibitor<-subset(query_result, moa =="Tyrosine kinase inhibitor")
query_result_Adenosine_deaminase_inhibitor<-subset(query_result, moa =="Adenosine deaminase inhibitor")
query_result_AKT_inhibitor<-subset(query_result, moa =="AKT inhibitor")
query_result_AMPK_inhibitor<-subset(query_result, moa =="AMPK inhibitor")
query_result_Angiogenesis_inhibitor<-subset(query_result, moa =="Angiogenesis inhibitor")
query_result_ATPase_inhibitor<-subset(query_result, moa =="ATPase inhibitor")
query_result_Cyclooxygenase_inhibitor<-subset(query_result, moa =="Cyclooxygenase inhibitor")
query_result_Focal_adhesion_kinase_inhibitor<-subset(query_result, moa =="Focal adhesion kinase inhibitor")
query_result_Guanylyl_cyclase_inhibitor<-subset(query_result, moa =="Guanylyl cyclase inhibitor")
query_result_HSP_inhibitor<-subset(query_result, moa =="HSP inhibitor")
query_result_IKK_inhibitor<-subset(query_result, moa =="IKK inhibitor")
query_result_JAK_inhibitor<-subset(query_result, moa =="JAK inhibitor")
query_result_Mediator_release_inhibitor<-subset(query_result, moa =="Mediator release inhibitor")
query_result_PARP_inhibitor<-subset(query_result, moa =="PARP inhibitor")
query_result_PDGFR_inhibitor<-subset(query_result, moa =="PDGFR inhibitor")
query_result_PKC_inhibitor<-subset(query_result, moa =="PKC inhibitor")
query_result_Protein_synthesis_inhibitor<-subset(query_result, moa =="Protein synthesis inhibitor")
query_result_Pyruvate_dehydrogenase_kinase_inhibitor<-subset(query_result, moa =="Pyruvate dehydrogenase kinase inhibitor")
query_result_Ribonucleotide_reductase_inhibitor<-subset(query_result, moa =="Ribonucleotide reductase inhibitor")
query_result_RNA_synthesis_inhibitor<-subset(query_result, moa =="RNA synthesis inhibitor")
query_result_Smoothened_receptor_antagonist<-subset(query_result, moa =="Smoothened receptor antagonist")

query_merged<-rbind(query_result_HDAC_inhibitor, query_result_Topoisomerase_inhibitor, query_result_CDK_inhibitor, 
      query_result_JNK_inhibitor, query_result_MTOR_inhibitor, query_result_RNA_polymerase_inhibitor, 
      query_result_Tyrosine_kinase_inhibitor, query_result_Adenosine_deaminase_inhibitor, query_result_AKT_inhibitor, 
      query_result_AMPK_inhibitor, query_result_Angiogenesis_inhibitor, query_result_ATPase_inhibitor, 
      query_result_Cyclooxygenase_inhibitor, query_result_Focal_adhesion_kinase_inhibitor, 
      query_result_Guanylyl_cyclase_inhibitor, query_result_HSP_inhibitor, query_result_IKK_inhibitor, 
      query_result_JAK_inhibitor, query_result_Mediator_release_inhibitor, query_result_PARP_inhibitor, 
      query_result_PDGFR_inhibitor, query_result_PKC_inhibitor, query_result_Protein_synthesis_inhibitor, 
      query_result_Pyruvate_dehydrogenase_kinase_inhibitor, query_result_Ribonucleotide_reductase_inhibitor, 
      query_result_RNA_synthesis_inhibitor, query_result_Smoothened_receptor_antagonist)

mydata<-subset(query_merged, select = c("pert_iname", "moa"))
write.csv(unique(mydata$moa),"moa.csv")
write.csv(unique(mydata$pert_iname),"pert_iname.csv")
t1<-unique(mydata$pert_iname) #构建行名
t2<-unique(mydata$moa) #构建列名
#214*27
datamatrix<-matrix(0,214,27)
rownames(datamatrix)<-t1
colnames(datamatrix)<-t2
for (i in 1:nrow(mydata)) {
  datamatrix[mydata[i,1],mydata[i,2]]<-1 #对应位置打上1
}
write.table(datamatrix,"datamatrix.txt",sep = "\t",quote = FALSE)
barplot_data<-colSums(datamatrix)
barplot_data<-as.data.frame(barplot_data)
barplot_data$moa<-row.names(barplot_data)
barplot_data$counts<-barplot_data$barplot_data
barplot_data<-subset(barplot_data, select = c("moa", "counts"))
barplot_data<-barplot_data[order(barplot_data$counts, decreasing = T),]
write.csv(barplot_data, "Barplot_dat.csv", row.names = F, quote = F)
pdf("barplot_MOAcounts.pdf")
barplot(barplot_data$counts)
dev.off()

merged_data<-merge(mydata, barplot_data, by = "moa")
merged_data_sort<-merged_data[order(merged_data$counts, decreasing = T),]

merged_data_sort$moa<-factor(merged_data_sort$moa, levels = barplot_data$moa)

library(ggplot2)
pdf("bubble_sorted.pdf", height = 20, width = 8)
p<-ggplot(merged_data_sort, aes(x = moa, y = pert_iname, size = 0.01)) + geom_point()+
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "gray"),
        panel.border = element_rect(colour="black",fill=NA))
p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = .5))
dev.off()

