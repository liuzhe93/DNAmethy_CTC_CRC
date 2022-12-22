setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/03_IndependentPrognosticAnalysis/02_Time-dependent_ROC")
rm(list=ls())

library(timeROC)
library(survival)
library(survivalROC)
###########################################TCGA数据的ROC曲线分析####################################################
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
merged_data$Gender_num<-ifelse(merged_data$Gender=="female","2",
                               ifelse(merged_data$Gender=="male","1",
                                      "NA"))
merged_data$Age_num<-ifelse(merged_data$Age>=60,"1",
                            ifelse(merged_data$Age<60,"2",
                                   "NA"))
merged_data$Stage_num<-ifelse(merged_data$Stage=="Stage IA", "1",
                              ifelse(merged_data$Stage=="Stage IIA", "2",
                                     ifelse(merged_data$Stage=="Stage IIB", "2",
                                            ifelse(merged_data$Stage=="Stage IIC", "2",
                                                   ifelse(merged_data$Stage=="Stage IIIA", "3",
                                                          ifelse(merged_data$Stage=="Stage IIIB", "3",
                                                                 ifelse(merged_data$Stage=="Stage IIIC", "3",
                                                                        ifelse(merged_data$Stage=="Stage IVA", "4",
                                                                               ifelse(merged_data$Stage=="Stage IVB", "4",
                                                                                      ifelse(merged_data$Stage=="Stage I", "1",
                                                                                             ifelse(merged_data$Stage=="Stage II", "2",
                                                                                                    ifelse(merged_data$Stage=="Stage III", "3",
                                                                                                           ifelse(merged_data$Stage=="Stage IV", "4",
                                                                                                                  "NA")))))))))))))

head(merged_data)
data_selected<-subset(merged_data, select = c("sample_name","futime","Status","TRIP10","NGFR","SLC48A1","SRMS","riskScore",
                                              "risk","Gender_num","Age_num","Stage_num"))
head(data_selected)
colnames(data_selected)<-c("sample_name","futime","Status","TRIP10","NGFR","SLC48A1","SRMS","riskScore","risk","Gender",
                           "Age","Stage")
head(data_selected)
write.csv(data_selected,"indepInput.csv",quote = F, row.names = F)
############################################ROC for RiskScore############################################################
time_roc_res_riskscore <- timeROC(
  T = data_selected$futime,
  delta = data_selected$Status,
  marker = data_selected$riskScore,
  cause = 1,
  weighting="marginal",
  times = c(1,2,3,4,5,6,7,8,9,10),
  ROC = TRUE,
  iid = TRUE
)
time_roc_res_riskscore$AUC
#t=1       t=2       t=3       t=4       t=5       t=6       t=7       t=8       t=9      t=10 
#0.7279297 0.7024452 0.6824516 0.6822810 0.6389583 0.5999614 0.6032997 0.5819867 0.6345913 0.6221901 
confint(time_roc_res_riskscore, level = 0.95)$CI_AUC
#      2.5% 97.5%
#t=1  66.06 79.53
#t=2  63.63 76.86
#t=3  60.52 75.97
#t=4  59.18 77.28
#t=5  53.56 74.23
#t=6  48.50 71.49
#t=7  47.88 72.78
#t=8  43.83 72.56
#t=9  47.87 79.05
#t=10 41.99 82.45
############################################ROC for Gender  ############################################################
data_selected$Gender<-as.numeric(data_selected$Gender)
time_roc_res_gender <- timeROC(
  T = data_selected$futime,
  delta = data_selected$Status,
  marker = data_selected$Gender,
  cause = 1,
  weighting="marginal",
  times = c(1,2,3,4,5,6,7,8,9,10),
  ROC = TRUE,
  iid = TRUE
)
time_roc_res_gender$AUC
#t=1       t=2       t=3       t=4       t=5       t=6       t=7       t=8       t=9      t=10 
#0.5023080 0.4811369 0.5007616 0.5000051 0.5101589 0.4725290 0.4446362 0.3957633 0.2770263 0.2603596 
confint(time_roc_res_gender, level = 0.95)$CI_AUC
#     2.5% 97.5%
#t=1  43.07 57.40
#t=2  41.71 54.52
#t=3  43.15 57.00
#t=4  41.55 58.45
#t=5  41.57 60.46
#t=6  36.33 58.18
#t=7  32.27 56.66
#t=8  24.88 54.27
#t=9  13.86 41.55
#t=10  9.90 42.17
############################################ROC for Age  ###############################################################
data_selected$Age<-as.numeric(data_selected$Age)
time_roc_res_age <- timeROC(
  T = data_selected$futime,
  delta = data_selected$Status,
  marker = data_selected$Age,
  cause = 1,
  weighting="marginal",
  times = c(1,2,3,4,5,6,7,8,9,10),
  ROC = TRUE,
  iid = TRUE
)
time_roc_res_age$AUC
#t=1       t=2       t=3       t=4       t=5       t=6       t=7       t=8       t=9      t=10 
#0.4606222 0.4905825 0.4922099 0.4169887 0.4017339 0.4610082 0.4796579 0.4547518 0.3944908 0.4278242 
confint(time_roc_res_age, level = 0.95)$CI_AUC
#      2.5% 97.5%
#t=1  40.23 51.89
#t=2  43.65 54.47
#t=3  43.48 54.96
#t=4  34.68 48.72
#t=5  31.79 48.56
#t=6  35.76 56.44
#t=7  37.05 58.88
#t=8  31.77 59.18
#t=9  23.44 55.46
#t=10 23.24 62.32
############################################ROC for Stage  #############################################################
data_selected$Stage<-as.numeric(data_selected$Stage)
time_roc_res_stage <- timeROC(
  T = data_selected$futime,
  delta = data_selected$Status,
  marker = data_selected$Stage,
  cause = 1,
  weighting="marginal",
  times = c(1,2,3,4,5,6,7,8,9,10),
  ROC = TRUE,
  iid = TRUE
)
time_roc_res_stage$AUC
#t=1       t=2       t=3       t=4       t=5       t=6       t=7       t=8       t=9      t=10 
#0.6861509 0.7015775 0.7297068 0.6972472 0.6517777 0.6718551 0.6300018 0.5918238 0.5106164 0.5106164
confint(time_roc_res_stage, level = 0.95)$CI_AUC
#     2.5% 97.5%
#t=1  61.08 76.15
#t=2  63.29 77.02
#t=3  66.23 79.71
#t=4  61.34 78.11
#t=5  55.32 75.04
#t=6  56.74 77.64
#t=7  51.80 74.20
#t=8  46.42 71.95
#t=9  36.59 65.54
#t=10 33.82 68.31
########################################绘图######################################################################
mydata<-cbind(time_roc_res_riskscore$AUC,time_roc_res_age$AUC,time_roc_res_gender$AUC,time_roc_res_stage$AUC)
colnames(mydata)<-c("Risk Score","Age", "Gender", "Stage")
mydata<-as.data.frame(mydata)
mydata$Year<-c(1,2,3,4,5,6,7,8,9,10)
all_data<-reshape2::melt(mydata, id.vars = "Year", variable.names = "DataSet")
write.csv(all_data, "AUC_timedependent.csv", row.names = F, quote = F)
#all_data<-read.csv("AUC_timedependent.csv",header=T)
library(ggplot2)
pdf("tROC_4factors.pdf")
theme_set(theme_bw())
p<-ggplot(data = all_data, mapping = aes(x = Year, y = value, colour = variable)) + geom_line() + ylim(c(0,1))
p_bottom=p+theme(legend.position = "top")
p_bottom
dev.off()

