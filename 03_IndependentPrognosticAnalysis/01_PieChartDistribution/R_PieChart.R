setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/03_IndependentPrognosticAnalysis/01_PieChartDistribution")
rm(list=ls())

#################################################TCGA数据分析##################################################################
risk<-read.table("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/04_multivariate_Cox_regression_analysis/risk.txt",
                 sep="\t",header=T, row.names = 1)
head(risk)
result_time<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/01_IdentificationDEGs/Metastatic_nonMetastatic/result_time.csv",header=T)
head(result_time)

risk$sample_name<-row.names(risk)
risk$sample_name<-substr(risk$sample_name,1,12)
head(risk)
merged_data<-merge(risk, result_time, by.x = "sample_name", by.y = "submitter_id")
merged_data<-subset(merged_data, select = c("sample_name","futime.x","fustat.x","TRIP10","NGFR","SLC48A1","SRMS","riskScore",
                                            "risk","gender","age_at_index","ajcc_pathologic_stage"))
colnames(merged_data)<-c("sample_name","futime","Status","TRIP10","NGFR","SLC48A1","SRMS","riskScore","risk","Gender",
                         "Age","Stage")
library("ggplot2")
library("ggforce")
merged_data_high<-subset(merged_data,risk=="high")
merged_data_low<-subset(merged_data,risk=="low")
###############################################TCGA_Status#################################################################
#变量：生存分析的变量有两个：生存时间t和结局变量(0-1)。其中结局变量1表示死亡事件，0表示截尾。
########high risk group##############
merged_data_high_status<-merged_data_high[,"Status"]
table(merged_data_high$Status)
#  0   1 
#187  72
num<-c(187,72)
type<-c("Alive","Dead")
status_high<-data.frame(num,type)
pdf("TCGA_high_status.pdf")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab("")+
  scale_fill_manual(values = c("red","green"))+
  geom_arc_bar(data = status_high, stat = "pie", aes(x0=0, y0=0, r0=1, r=2,amount=num, fill=type))+
  theme(legend.position = "top")
dev.off()
########low risk group##############
merged_data_low_status<-merged_data_low[,"Status"]
table(merged_data_low$Status)
#  0   1 
#217  43
num<-c(217,43)
type<-c("Alive","Dead")
status_low<-data.frame(num,type)
pdf("TCGA_low_status.pdf")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab("")+
  scale_fill_manual(values = c("red","green"))+
  geom_arc_bar(data = status_low, stat = "pie", aes(x0=0, y0=0, r0=1, r=2,amount=num, fill=type))+
  theme(legend.position = "top")
dev.off()
merged_status_data<-cbind(status_high,status_low)
row.names(merged_status_data)<-merged_status_data$type
merged_status_data<-merged_status_data[,-2]
merged_status_data<-merged_status_data[,-3]
colnames(merged_status_data)<-c("high", "low")
chisq.test(merged_status_data)
#Pearson's Chi-squared test with Yates' continuity correction
#
#data:  merged_status_data
#X-squared = 8.8972, df = 1, p-value = 0.002856


###############################################TCGA_Stage#################################################################
merged_data_cp<-merged_data
merged_data$Stage_ope<-ifelse(merged_data$Stage=="Stage IA", "Stage I",
                              ifelse(merged_data$Stage=="Stage IIA", "Stage II",
                                     ifelse(merged_data$Stage=="Stage IIB", "Stage II",
                                            ifelse(merged_data$Stage=="Stage IIC", "Stage II",
                                                   ifelse(merged_data$Stage=="Stage IIIA", "Stage III",
                                                          ifelse(merged_data$Stage=="Stage IIIB", "Stage III",
                                                                 ifelse(merged_data$Stage=="Stage IIIC", "Stage III",
                                                                        ifelse(merged_data$Stage=="Stage IVA", "Stage IV",
                                                                               ifelse(merged_data$Stage=="Stage IVB", "Stage IV",
                                                                                      ifelse(merged_data$Stage=="Stage I", "Stage I",
                                                                                             ifelse(merged_data$Stage=="Stage II", "Stage II",
                                                                                                    ifelse(merged_data$Stage=="Stage III", "Stage III",
                                                                                                           ifelse(merged_data$Stage=="Stage IV", "Stage IV",
                                                                                                                  "NA")))))))))))))
merged_data_high<-subset(merged_data,risk=="high")
merged_data_low<-subset(merged_data,risk=="low")
########high risk group##############
merged_data_high_stage<-merged_data_high[,"Stage_ope"]
table(merged_data_high$Stage_ope)
#Stage I  Stage II Stage III  Stage IV 
#     34        93        76        51
num<-c(34,93,76,51)
type<-c("Stage I","Stage II","Stage III","Stage IV")
stage_high<-data.frame(num,type)
pdf("TCGA_high_stage.pdf")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab("")+
  scale_fill_manual(values = c("red","green","blue","black"))+
  geom_arc_bar(data = stage_high, stat = "pie", aes(x0=0, y0=0, r0=1, r=2,amount=num, fill=type))+
  theme(legend.position = "top")
dev.off()
########low risk group##############
merged_data_low_stage<-merged_data_low[,"Stage_ope"]
table(merged_data_low$Stage_ope)
#Stage I  Stage II Stage III  Stage IV 
#     51       116        64        22
num<-c(51,116,64,22)
type<-c("Stage I","Stage II","Stage III","Stage IV")
stage_low<-data.frame(num,type)
pdf("TCGA_low_stage.pdf")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab("")+
  scale_fill_manual(values = c("red","green","blue","black"))+
  geom_arc_bar(data = stage_low, stat = "pie", aes(x0=0, y0=0, r0=1, r=2,amount=num, fill=type))+
  theme(legend.position = "top")
dev.off()
merged_stage_data<-cbind(stage_high,stage_low)
row.names(merged_stage_data)<-merged_stage_data$type
merged_stage_data<-merged_stage_data[,-2]
merged_stage_data<-merged_stage_data[,-3]
colnames(merged_stage_data)<-c("high", "low")
chisq.test(merged_stage_data)
#Pearson's Chi-squared test
#
#data:  merged_stage_data
#X-squared = 18.478, df = 3, p-value = 0.0003504
###############################################TCGA_Gender#################################################################
#变量：生存分析的变量有两个：生存时间t和结局变量(0-1)。其中结局变量1表示死亡事件，0表示截尾。
########high risk group##############
merged_data_high_Gender<-merged_data_high[,"Gender"]
table(merged_data_high$Gender)
#female   male 
#125    134
num<-c(125,134)
type<-c("Female","Male")
gender_high<-data.frame(num,type)
pdf("TCGA_high_gender.pdf")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab("")+
  scale_fill_manual(values = c("red","green"))+
  geom_arc_bar(data = gender_high, stat = "pie", aes(x0=0, y0=0, r0=1, r=2,amount=num, fill=type))+
  theme(legend.position = "top")
dev.off()
########low risk group##############
merged_data_low_Gender<-merged_data_low[,"Gender"]
table(merged_data_low$Gender)
#female   male 
#122    138
num<-c(122,138)
type<-c("Female","Male")
gender_low<-data.frame(num,type)
pdf("TCGA_low_gender.pdf")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab("")+
  scale_fill_manual(values = c("red","green"))+
  geom_arc_bar(data = gender_low, stat = "pie", aes(x0=0, y0=0, r0=1, r=2,amount=num, fill=type))+
  theme(legend.position = "top")
dev.off()
merged_gender_data<-cbind(gender_high,gender_low)
row.names(merged_gender_data)<-merged_gender_data$type
merged_gender_data<-merged_gender_data[,-2]
merged_gender_data<-merged_gender_data[,-3]
colnames(merged_gender_data)<-c("high", "low")
chisq.test(merged_gender_data)
#Pearson's Chi-squared test with Yates' continuity correction
#
#data:  merged_gender_data
#X-squared = 0.047356, df = 1, p-value = 0.8277

###############################################TCGA_Age#################################################################
#变量：生存分析的变量有两个：生存时间t和结局变量(0-1)。其中结局变量1表示死亡事件，0表示截尾。
merged_data$Age_ope<-ifelse(merged_data$Age>=60,">=60","<60")
merged_data_high<-subset(merged_data,risk=="high")
merged_data_low<-subset(merged_data,risk=="low")
########high risk group##############
merged_data_high_age<-merged_data_high[,"Age_ope"]
table(merged_data_high$Age_ope)
#<60 >=60 
#73  186 
num<-c(73,186)
type<-c("<60",">=60")
age_high<-data.frame(num,type)
pdf("TCGA_high_age.pdf")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab("")+
  scale_fill_manual(values = c("red","green"))+
  geom_arc_bar(data = age_high, stat = "pie", aes(x0=0, y0=0, r0=1, r=2,amount=num, fill=type))+
  theme(legend.position = "top")
dev.off()
########low risk group##############
merged_data_low_age<-merged_data_low[,"Age_ope"]
table(merged_data_low$Age_ope)
#<60 >=60 
#63  197
num<-c(63,197)
type<-c("<60",">=60")
age_low<-data.frame(num,type)
pdf("TCGA_low_age.pdf")
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab("")+
  scale_fill_manual(values = c("red","green"))+
  geom_arc_bar(data = age_low, stat = "pie", aes(x0=0, y0=0, r0=1, r=2,amount=num, fill=type))+
  theme(legend.position = "top")
dev.off()
merged_age_data<-cbind(age_high,age_low)
row.names(merged_age_data)<-merged_age_data$type
merged_age_data<-merged_age_data[,-2]
merged_age_data<-merged_age_data[,-3]
colnames(merged_age_data)<-c("high", "low")
chisq.test(merged_age_data)
#Pearson's Chi-squared test with Yates' continuity correction
#
#data:  merged_age_data
#X-squared = 0.85476, df = 1, p-value = 0.3552


