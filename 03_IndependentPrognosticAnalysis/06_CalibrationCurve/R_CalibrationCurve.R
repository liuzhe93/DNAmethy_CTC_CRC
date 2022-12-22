setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/03_IndependentPrognosticAnalysis/06_CalibrationCurve")
rm(list = ls())

library(rms)
library(timeROC)
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
                                              "risk","Gender","Age","Stage_num"))
colnames(data_selected)<-c("sample_name","futime","Status","TRIP10","NGFR","SLC48A1","SRMS","riskScore","risk","Gender",
                           "Age","Stage")
rt<-data_selected
head(rt)
for(i in 1:(dim(rt)[1])){
  rt$sample_name_new[i]<-paste(rt$sample_name[i],i,sep="_")
}
rownames(rt)<-rt$sample_name_new
rt<-subset(rt, select = c("futime", "Status", "Age", "Gender", "Stage", "riskScore"))
colnames(rt)<-c("futime", "fustat", "Age", "Gender", "Stage", "riskScore")
rt$Stage<-ifelse(rt$Stage=="1","Stage I",
                 ifelse(rt$Stage=="2","Stage II",
                        ifelse(rt$Stage=="3","Stage III",
                               ifelse(rt$Stage=="4","Stage IV",
                                      "NA"))))
head(rt)
ddist <- datadist(rt)
options(datadist='ddist')
f_cph <- cph(Surv(futime,fustat) ~ Age+Gender+Stage+riskScore,
             x=T, y=T, surv=T,
             data=rt)
#查看多因素Cox分析结果，最下方可见其对应的Coef、p值
print(f_cph)
ddist <- datadist(rt)
options(datadist='ddist')
med  <- Quantile(f_cph)
surv <- Survival(f_cph) 
cal_1<-calibrate(f_cph,u=3,cmethod='KM',m=15,B=200)# usually B=200 or 300
cal_2<-calibrate(f_cph,u=5,cmethod='KM',m=15,B=200)# usually B=200 or 300
cal_3<-calibrate(f_cph,u=10,cmethod='KM',m=15,B=200)# usually B=200 or 300
#par(mar=c(7,4,4,3),cex=1.0)
pdf("calibrate_3_years.pdf",width=6,height=6) 
plot(cal_1,lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1 years OS',#便签
     ylab='Actual 3 years OS(proportion)',#标签
     col=c(rgb(192,98,83,maxColorValue = 255)),#设置一个颜色
     xlim = c(0,1),ylim = c(0,1),##x轴和y轴范围
     mgp = c(2, 1, 0)) #控制坐标轴的位置
dev.off()
pdf("calibrate_5_years.pdf",width=6,height=6) 
plot(cal_2,lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 5 years OS',#便签
     ylab='Actual 5 years OS(proportion)',#标签
     col=c(rgb(192,98,83,maxColorValue = 255)),#设置一个颜色
     xlim = c(0,1),ylim = c(0,1),##x轴和y轴范围
     mgp = c(2, 1, 0)) #控制坐标轴的位置
dev.off()
pdf("calibrate_10_years.pdf",width=6,height=6) 
plot(cal_3,lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 10 years OS',#便签
     ylab='Actual 10 years OS(proportion)',#标签
     col=c(rgb(192,98,83,maxColorValue = 255)),#设置一个颜色
     xlim = c(0,1),ylim = c(0,1),##x轴和y轴范围
     mgp = c(2, 1, 0)) #控制坐标轴的位置
dev.off()

