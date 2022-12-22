setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/07_ChemotherapyResistanceAnalysis/GSE143985/drugs")
rm(list = ls())

risk<-read.table("D:/MyProjects/CTC/CTC_Prognosis_COAD/07_ChemotherapyResistanceAnalysis/GSE143985/risk.txt",
                 header=T,row.names = 1)
risk$sample_name<-row.names(risk)
dim(risk)
#[1] 91  9
ic50<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/07_ChemotherapyResistanceAnalysis/GSE143985_drugSensitivity.csv")
ic50$sample_name<-ic50$X
dim(ic50)
#[1] 91  62

merged<-merge(ic50,risk,by="sample_name")
dim(merged)
#[1] 91  70

rownames(merged)<-merged$X
merged_format<-cbind(merged[,ncol(merged)],merged[,3:(ncol(merged)-2)])
####################################drug1： Cisplatin################################################
drug_sen_Cisplatin<-cbind(merged_format$`merged[, ncol(merged)]`,as.numeric(merged_format$Cisplatin))
rownames(drug_sen_Cisplatin)<-row.names(merged_format)
colnames(drug_sen_Cisplatin)<-c("Risk","IC50")
head(drug_sen_Cisplatin)
drug_sen_Cisplatin<-as.data.frame(drug_sen_Cisplatin)
drug_sen_Cisplatin$IC50<-as.numeric(drug_sen_Cisplatin$IC50)
library(ggpubr)
pdf("Cisplatin.pdf")
p <- ggboxplot(drug_sen_Cisplatin, x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug2： Camptothecin ################################################
drug_sen_Camptothecin <-cbind(merged_format$`merged[, ncol(merged)]`,as.numeric(merged_format$Camptothecin))
rownames(drug_sen_Camptothecin)<-row.names(merged_format)
colnames(drug_sen_Camptothecin)<-c("Risk","IC50")
head(drug_sen_Camptothecin)
drug_sen_Camptothecin <-as.data.frame(drug_sen_Camptothecin)
drug_sen_Camptothecin$IC50<-as.numeric(drug_sen_Camptothecin$IC50)
pdf("Camptothecin.pdf")
p <- ggboxplot(drug_sen_Camptothecin , x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()

####################################drug3： Rapamycin ################################################
drug_sen_Rapamycin <-cbind(merged_format$`merged[, ncol(merged)]`,as.numeric(merged_format$Rapamycin))
rownames(drug_sen_Rapamycin )<-row.names(merged_format)
colnames(drug_sen_Rapamycin )<-c("Risk","IC50")
head(drug_sen_Rapamycin )
drug_sen_Rapamycin <-as.data.frame(drug_sen_Rapamycin )
drug_sen_Rapamycin $IC50<-as.numeric(drug_sen_Rapamycin $IC50)
pdf("Rapamycin .pdf")
p <- ggboxplot(drug_sen_Rapamycin , x="Risk",
               y = "IC50", color = "Risk",
               palette = "jco", add = "jitter")
# 添加p值
p + stat_compare_means()
dev.off()











library(ggpubr)
library(pRRophetic)
library(ggplot2)
data(cgp2016ExprRma) 
dim(cgp2016ExprRma)
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)
possibleDrugs2016 <- unique(drugData2016$Drug.name)
remove_drug<-c("PHA-665752","MG-132","VX-680","TAE684","Crizotinib","Saracatinib","S-Trityl-L-cysteine",
               "Z-LLNle-CHO","GNF-2","CGP-60474","CGP-082996","A-770041","WH-4-023","WZ-1-84","BI-2536",
               "BMS-536924","BMS-509744","JW-7-52-1","A-443654","MS-275","KIN001-135","TGX221","XMD8-85",
               "Mitomycin C","NSC-87877","CP466722","CHIR-99021","AP-24534","JNK-9L","PF-562271","HG-6-64-1",
               "JQ1","JQ12","FTI-277","OSU-03012","AKT inhibitor VIII","PAC-1","IPA-3","GSK-650394",
               "BAY 61-3606","5-Fluorouracil","Obatoclax Mesylate","BMS-754807","Lisitinib","LFM-A13",
               "GW-2580","Phenformin","Bryostatin 1","LAQ824","Epothilone B","GSK1904529A","BMS345541",
               "BMS-708163","Ruxolitinib","Ispinesib Mesylate","TL-2-105","AT-7519","TAK-715","BX-912",
               "ZSTK474","AS605240","Genentech Cpd 10","GSK1070916","KIN001-102","LY317615","GSK429286A","FMK",
               "QL-XII-47","CAL-101","UNC0638","XL-184","WZ3105","XMD14-99","AC220","CP724714","JW-7-24-1",
               "NPK76-II-72-1","STF-62247","NG-25","TL-1-85","VX-11e","FR-180204","Tubastatin A","Zibotentan",
               "YM155","NSC-207895","VNLG/124","AR-42","CUDC-101","Belinostat","I-BET-762","CAY10603",
               "Linifanib ","BIX02189","CH5424802","EKB-569","GSK2126458","KIN001-236","KIN001-244",
               "KIN001-055","KIN001-260","KIN001-266","Masitinib","MP470","MPS-1-IN-1","BHG712","OSI-930",
               "OSI-027","CX-5461","PHA-793887","PI-103","PIK-93","SB52334","TPCA-1","TG101348","Foretinib",
               "Y-39983","YM201636","Tivozanib","GSK690693","SNX-2112","QL-XI-92","XMD13-2","QL-X-138",
               "XMD15-27","T0901317","EX-527","THZ-2-49","KIN001-270","THZ-2-102-1","Navitoclax","CI-1040",
               "Olaparib","Veliparib","VX-702","AMG-706","KU-55933","Afatinib","GDC0449","BX-795","NU-7441",
               "SL 0101-1","BIRB 0796","JNK Inhibitor VIII","681640","Nutlin-3a (-)","PD-173074","ZM-447439",
               "RO-3306","MK-2206","PD-0332991","BEZ235","PD-0325901","selumetinib","EHT 1864","Cetuximab",
               "PF-4708671","JNJ-26854165","HG-5-113-01","HG-5-88-01","TW 37","XMD11-85h","ZG-10","XMD8-92",
               "QL-VIII-58","AG-014699","SB 505124","Tamoxifen","QL-XII-61","PFI-1","IOX2","YK 4-279",
               "(5Z)-7-Oxozeaenol","piperlongumine","FK866","Talazoparib","rTRAIL","UNC1215","SGC0946",
               "XAV939","Trametinib","Dabrafenib","Temozolomide","Bleomycin (50 uM)","SN-38","MLN4924",
               "AZD7762","GW 441756","CEP-701","SB 216763","17-AAG")
drug_list<-setdiff(possibleDrugs2016,remove_drug)
write.csv(drug_list, "Drugs.csv", quote = F)
setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/07_ChemotherapyResistanceAnalysis/GSE143985/drugs/Figures/")
for(drug_each in drug_list){
  print(drug_each)
  drug_sen<-cbind(merged_format$`merged[, ncol(merged)]`, subset(merged_format, select = drug_each))
  rownames(drug_sen)<-row.names(merged_format)
  colnames(drug_sen)<-c("Risk","IC50")
  head(drug_sen)
  drug_sen<-as.data.frame(drug_sen)
  drug_sen$IC50<-as.numeric(drug_sen$IC50)
  varnames<-paste(drug_each,".pdf",sep = "")
  pdf(varnames)
  p <- ggboxplot(drug_sen , x="Risk",
                 y = "IC50", color = "Risk",
                 palette = "jco", add = "jitter")
  # 添加p值
  p + stat_compare_means()
  dev.off()
  
}
