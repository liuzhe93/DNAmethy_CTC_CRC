setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/06_GenomicVariation/SNV")
rm(list = ls()) 

library("TCGAbiolinks")
library("tidyverse")
library("maftools")

####################################read clinical data#####################################################
clin_data<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/Metastatic_nonMetastatic/result_time.csv",header=T)
head(clin_data)

####################################risk Score data########################################################
risk<-read.table("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/04_multivariate_Cox_regression_analysis/risk.txt",
                 sep="\t",header=T, row.names = 1)
head(risk)
risk$sample_name<-row.names(risk)
risk$sample_name<-substr(risk$sample_name,1,12)
head(risk)
merged_data<-merge(risk, clin_data, by.x = "sample_name", by.y = "submitter_id")
merged_data<-subset(merged_data, select = c("sample_name","futime.x","fustat.x","TRIP10","NGFR","SLC48A1","SRMS","riskScore",
                                            "risk","gender","age_at_index","ajcc_pathologic_stage"))
colnames(merged_data)<-c("sample_name","futime","Status","TRIP10","NGFR","SLC48A1","SRMS","riskScore","risk","Gender",
                         "Age","Stage")
library("ggplot2")
library("ggforce")
merged_data_high<-subset(merged_data,risk=="high")
merged_data_low<-subset(merged_data,risk=="low")

####################################read SNP data##########################################################
maf.data<-read.table("SNP/maf.tsv", header = T)
selcol=c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type" , "Tumor_Sample_Barcode")
maftools_df=maf.data[,selcol]
write.table(data.frame(maftools_df,check.names = F), file = 'COAD_maftools_df.maf', sep="\t", row.names =F, quote = F)
maf = read.maf(maf ='COAD_maftools_df.maf')
maftools_df$sample_name<-maftools_df$Tumor_Sample_Barcode
maftools_df$sample_name<-substr(maftools_df$sample_name,1,12)
merged_maf<-merge(maftools_df,risk,by="sample_name")
merged_maf_high<-subset(merged_maf, risk == "high")
merged_maf_high<-subset(merged_maf_high, select = c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", 
                                               "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type" , "Tumor_Sample_Barcode"))
write.table(data.frame(merged_maf_high,check.names = F), file = 'COAD_maftools_df_high.maf', sep="\t", row.names =F, quote = F)
maf_high = read.maf(maf ='COAD_maftools_df_high.maf')

merged_maf_low<-subset(merged_maf, risk == "low")
merged_maf_low<-subset(merged_maf_low, select = c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", 
                                                    "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type" , "Tumor_Sample_Barcode"))
write.table(data.frame(merged_maf_low,check.names = F), file = 'COAD_maftools_df_low.maf', sep="\t", row.names =F, quote = F)
maf_low = read.maf(maf ='COAD_maftools_df_low.maf')

pdf("WaterfallPlots_high.pdf")
oncoplot(maf = maf_high, top = 30)
dev.off()
pdf("WaterfallPlots_low.pdf")
oncoplot(maf = maf_low, top = 30)
dev.off()

clindata = maf@clinical.data ###提取SNP中对应病人的临床数据
plotmafSummary(maf, rmOutlier = TRUE, dashboard = TRUE,
               titvRaw = TRUE, addStat = NULL, showBarcodes = FALSE, fs = 0.8,
               textSize = 0.8, color = NULL, titleSize = c(0.8, 0.6),
               titvColor = NULL, top = 10)

Clin = getClinicalData(maf)
colNM = getFields(maf)
mafs <- mafSummary(maf)
lamaGene = getGeneSummary(maf)
lamaSample = getSampleSummary(maf)

############################################################################################################
lamaSample<-as.data.frame(lamaSample)
lamaSample$sample_name<-lamaSample$Tumor_Sample_Barcode
lamaSample$sample_name<-substr(lamaSample$sample_name, 1, 12)
lamaSample_merged<-merge(lamaSample,merged_data,by="sample_name")
lamaSample_merged_high<-subset(lamaSample_merged, risk == "high")
lamaSample_merged_low<-subset(lamaSample_merged, risk == "low")
dim(lamaSample_merged_high)
#[1] 254  23
dim(lamaSample_merged_low)
#[1] 292  23
for (i in 1:nrow(lamaSample_merged)) {
  lamaSample_merged$sample_name[i]<-paste(lamaSample_merged$sample_name[i],i,sep="_")
}

############################total mutation counts###########################################################
all_mutation<-subset(lamaSample_merged, select = c("sample_name","total","riskScore","risk"))
rownames(all_mutation)<-all_mutation$sample_name
all_mutation<-all_mutation[,-1]
all_mutation$total<-log2(all_mutation$total)
all_mutation_high<-subset(all_mutation, risk == "high")
all_mutation_low<-subset(all_mutation, risk == "low")
pdf("total.pdf")
p<-ggscatter(all_mutation,x="riskScore",y="total",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
pdf("total_high.pdf")
p<-ggscatter(all_mutation_high,x="riskScore",y="total",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
pdf("total_low.pdf")
p<-ggscatter(all_mutation_low,x="riskScore",y="total",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
################################Missense_Mutation###########################################################
mis_mutation<-subset(lamaSample_merged, select = c("sample_name","Missense_Mutation","riskScore","risk"))
rownames(mis_mutation)<-mis_mutation$sample_name
mis_mutation<-mis_mutation[,-1]
mis_mutation$Missense_Mutation<-log2(mis_mutation$Missense_Mutation)
mis_mutation_high<-subset(mis_mutation, risk == "high")
mis_mutation_low<-subset(mis_mutation, risk == "low")
pdf("Missense_Mutation.pdf")
p<-ggscatter(mis_mutation,x="riskScore",y="Missense_Mutation",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
pdf("Missense_Mutation_high.pdf")
p<-ggscatter(mis_mutation_high,x="riskScore",y="Missense_Mutation",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
pdf("Missense_Mutation_low.pdf")
p<-ggscatter(mis_mutation_low,x="riskScore",y="Missense_Mutation",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()

################################Nonsense_Mutation###########################################################
non_mutation<-subset(lamaSample_merged, select = c("sample_name","Nonsense_Mutation","riskScore","risk"))
rownames(non_mutation)<-non_mutation$sample_name
non_mutation<-non_mutation[,-1]
non_mutation$Nonsense_Mutation<-log2(non_mutation$Nonsense_Mutation)
non_mutation_high<-subset(non_mutation, risk == "high")
non_mutation_low<-subset(non_mutation, risk == "low")
pdf("Nonsense_Mutation.pdf")
p<-ggscatter(non_mutation,x="riskScore",y="Nonsense_Mutation",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
pdf("Nonsense_Mutation_high.pdf")
p<-ggscatter(non_mutation_high,x="riskScore",y="Nonsense_Mutation",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()
pdf("Nonsense_Mutation_low.pdf")
p<-ggscatter(non_mutation_low,x="riskScore",y="Nonsense_Mutation",
             add = "reg.line", conf.int = T,
             add.params = list(fill = "lightgray"))+
  stat_cor(method = "pearson")
p1<-ggMarginal(p+ggtitle("type='density'"), type = "density",
               xparams = list(fill ="orange"),
               yparams = list(fill ="skyblue"))
p1
dev.off()

#############################################
genelist_high <- c("APC", "TP53", "TTN", "KRAS", "SYNE1", "PIK3CA", "MUC16", "FAT4", "CSMD3", "DNAH5", 
                  "RYR2", "ZFHX4", "OBSCN", "USH2A", "RYR3", "LRP1B", "NEB", "RYR1", "ABCA13", "ADGRV1", 
                  "FBXW7", "KMT2D", "PCLO", "RNF43", "BRAF", "CSMD1", "FLG", "FAT3", "DST", "SYNE2")
genelist_low <- c("APC", "TP53", "TTN", "KRAS", "MUC16", "SYNE1", "PIK3CA", "FAT4", "RYR2", "OBSCN", 
                  "CSMD3", "DNAH11", "DNAH5", "DOCK3", "LRP1B", "ZFHX4", "FAT3", "PCLO", "CSMD1", "ABCA13", 
                  "ATM", "DST", "RYR3", "ANK3", "NEB", "ACVR2A", "CACNA1E", "DMD", "MDN1", "MUC5B")
overlapping_genes<-intersect(genelist_high, genelist_low)
lamaGene_high = getGeneSummary(maf_high)
lamaGene_high<-as.data.frame(lamaGene_high)
row.names(lamaGene_high)<-lamaGene_high$Hugo_Symbol
lamaGene_high_selected<-lamaGene_high[overlapping_genes,]
lamaGene_high_selected<-subset(lamaGene_high_selected, select = c("MutatedSamples"))

lamaGene_low = getGeneSummary(maf_low)
lamaGene_low<-as.data.frame(lamaGene_low)
row.names(lamaGene_low)<-lamaGene_low$Hugo_Symbol
lamaGene_low_selected<-lamaGene_low[overlapping_genes,]
lamaGene_low_selected<-subset(lamaGene_low_selected, select = c("MutatedSamples"))

##########################example#######################################################################################
library("fmsb")
dat<-matrix(c(23,6,117,210), nrow = 2, ncol = 2)
rownames(dat)<-c("Mutated genes", "Non-mutated genes")
colnames(dat)<-c("Cancer", "Normal")
fmsb::oddsratio(dat)
###############################Forest plot##############################################################################
mut_dat<-cbind(lamaGene_high_selected, lamaGene_low_selected)
colnames(mut_dat)<-c("high_mut_samples", "low_mut_samples")
write.csv(mut_dat, "ForestData_OR_Pvalue.csv", quote = F)
outTab=data.frame()
for(i in 1:nrow(mut_dat)){
  dat<-matrix(c(mut_dat[i,1],(254-mut_dat[i,1]),mut_dat[i,2],(292-mut_dat[i,2])), nrow = 2, ncol = 2)
  rownames(dat)<-c("Mutated genes", "Non-mutated genes")
  colnames(dat)<-c("HRisk", "LRisk")
  res<-fmsb::oddsratio(dat)
  outTab=rbind(outTab,
               cbind(gene=row.names(mut_dat)[i],
                     HR=res$estimate,
                     HR.95L=res$conf.int[1],
                     HR.95H=res$conf.int[2],
                     pvalue=res$p.value)
  )
}
write.csv(outTab,"HR_data.csv",quote=F,row.names = F)
rt<-read.csv("HR_data.csv", row.names = 1, header = T)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
#HRisk=mut_dat[,"high_mut_samples"]
#LRisk=mut_dat[,"low_mut_samples"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
#       append("LRisk", LRisk),
#       append("HRisk", HRisk),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )
library(forestplot)
pdf(file="forest.pdf",onefile = FALSE,
    width = 6,             #图片的宽度
    height = 4,            #图片的高度
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()
#######################################colinearity of mutations analysis##############################################
#download data from cbioportal via the following url
#https://www.cbioportal.org/results/mutualExclusivity?cancer_study_list=coadread_tcga&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&profileFilter=mutations%2Cgistic&case_set_id=coadread_tcga_cnaseq&gene_list=APC%250ATP53%250ATTN%250AKRAS%250ASYNE1%250APIK3CA%250AMUC16%250AFAT4%250ACSMD3%250ADNAH5%250ARYR2%250AZFHX4%250AOBSCN%250ARYR3%250ALRP1B%250ANEB%250AABCA13%250APCLO%250ACSMD1%250AFAT3%250ADST&geneset_list=%20&tab_index=tab_visualize&Action=Submit
overlapping_genes
# [1] "APC"    "TP53"   "TTN"    "KRAS"   "SYNE1"  "PIK3CA" "MUC16"  "FAT4"   "CSMD3"  "DNAH5"  "RYR2"   "ZFHX4"  "OBSCN"  "RYR3"   "LRP1B" 
#[16] "NEB"    "ABCA13" "PCLO"   "CSMD1"  "FAT3"   "DST"
write.csv(overlapping_genes, "overlapp21genes.csv", quote = F)
mydata_bio<-read.table("table.tsv", header = T, sep = "\t")
mydata_APC<-subset(mydata_bio, A == "APC")
mydata_TP53<-subset(mydata_bio, A == "TP53")
mydata_TTN<-subset(mydata_bio, A == "TTN")
mydata_KRAS<-subset(mydata_bio, A == "KRAS")
mydata_SYNE1<-subset(mydata_bio, A == "SYNE1")
mydata_PIK3CA<-subset(mydata_bio, A == "PIK3CA")
mydata_MUC16<-subset(mydata_bio, A == "MUC16")
mydata_FAT4<-subset(mydata_bio, A == "FAT4")
mydata_CSMD3<-subset(mydata_bio, A == "CSMD3")
mydata_DNAH5<-subset(mydata_bio, A == "DNAH5")
mydata_RYR2<-subset(mydata_bio, A == "RYR2")
mydata_ZFHX4<-subset(mydata_bio, A == "ZFHX4")
mydata_OBSCN<-subset(mydata_bio, A == "OBSCN")
mydata_RYR3<-subset(mydata_bio, A == "RYR3")
mydata_LRP1B<-subset(mydata_bio, A == "LRP1B")
mydata_NEB<-subset(mydata_bio, A == "NEB")
mydata_ABCA13<-subset(mydata_bio, A == "ABCA13")
mydata_PCLO<-subset(mydata_bio, A == "PCLO")
mydata_CSMD1<-subset(mydata_bio, A == "CSMD1")
mydata_FAT3<-subset(mydata_bio, A == "FAT3")

colinear_p<-read.csv("Cor_overlapp21genes.csv", header = T, row.names = 1)
library(corrplot)
colinear_p[colinear_p == "<0.001"]<-"0.0005"
colinear_p$APC<-as.numeric(colinear_p$APC)
colinear_p$TP53<-as.numeric(colinear_p$TP53)
colinear_p$TTN<-as.numeric(colinear_p$TTN)
colinear_p$KRAS<-as.numeric(colinear_p$KRAS)
colinear_p$SYNE1<-as.numeric(colinear_p$SYNE1)
colinear_p$PIK3CA<-as.numeric(colinear_p$PIK3CA)
colinear_p$MUC16<-as.numeric(colinear_p$MUC16)
colinear_p$FAT4<-as.numeric(colinear_p$FAT4)
colinear_p$CSMD3<-as.numeric(colinear_p$CSMD3)
colinear_p$DNAH5<-as.numeric(colinear_p$DNAH5)
colinear_p$RYR2<-as.numeric(colinear_p$RYR2)
colinear_p$ZFHX4<-as.numeric(colinear_p$ZFHX4)
colinear_p$OBSCN<-as.numeric(colinear_p$OBSCN)
colinear_p$RYR3<-as.numeric(colinear_p$RYR3)
colinear_p$LRP1B<-as.numeric(colinear_p$LRP1B)
colinear_p$NEB<-as.numeric(colinear_p$NEB)
colinear_p$ABCA13<-as.numeric(colinear_p$ABCA13)
colinear_p$PCLO<-as.numeric(colinear_p$PCLO)
colinear_p$CSMD1<-as.numeric(colinear_p$CSMD1)
colinear_p$FAT3<-as.numeric(colinear_p$FAT3)
colinear_p$DST<-as.numeric(colinear_p$DST)
colinear_p<-as.matrix(colinear_p)

library("pheatmap")
pdf("Heatmap_colinear_Legend.pdf", width = 8, height = 6)
pheatmap(-log2(colinear_p+0.001), cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
dev.off()

