setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/04_FunctionalEnrichmentAnalysis/03_CNVanalysis")
rm(list = ls())

# https://www.cbioportal.org/study/summary?id=coadread_tcga
# https://www.cbioportal.org/
# Liver(left of the page)-->  bowel --> Colorectal Adenocarcinoma (TCGA, Firehouse Legancy)
# 640 samples
# download CNA_Genes.txt file
cna_data<-read.table("CNA_Genes.txt",sep="\t",header = T)
geneList <- c("TRIP10", "NGFR", "SLC48A1", "SRMS")


gene_index<-t(t(geneList))
colnames(gene_index)<-"Gene"
gene_index<-as.data.frame(gene_index)
data_export<-merge(cna_data, gene_index, by = "Gene")
amp_freq<-subset(data_export,CNA=="AMP")
summary(as.numeric(sub("%","",amp_freq$Freq))/100)
del_freq<-subset(data_export,CNA=="HOMDEL")
summary(as.numeric(sub("%","",del_freq$Freq))/100)
del_freq$Freq<-paste0("-",del_freq$Freq)
amp_freq$Freq<-as.numeric(sub("%","",amp_freq$Freq))
del_freq$Freq<-as.numeric(sub("%","",del_freq$Freq))
data<-rbind(amp_freq,del_freq)

pdf("CNA.pdf",width=8,height = 6)
ggplot(data = data) +
  geom_col(aes(x = factor(Gene), y = Freq, fill = CNA)) +
  scale_y_continuous(breaks = seq(from = -8, to = 16,by = 1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = .5))
dev.off()

