setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/08_ImmunotherapyResponse")
rm(list=ls())
library(nnet)
data<-read.table("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/04_multivariate_Cox_regression_analysis/risk.txt",
                 sep = "\t", header = T)
#data$group<-ifelse(data$risk=="low",0,1)
#cha=0
#wine=data
#for(i in 1:521)
#{ 
#  if(wine[i,10]==0)cha[i]="low"
#  else cha[i]="high"
#}
#wine[,11]=factor(cha)
data$risk<-factor(data$risk)
scale01=function(x)
{
  ncol=dim(x)[2]-1#提取特征变量个数
  nrow=dim(x)[1]#提取样本集中的样本总量
  new=matrix(0,nrow,ncol)#建立用于保存新样本集的矩阵
  for(i in 1:ncol)
  {
    max=max(x[,i])#提取每个变量的最大值
    min=min(x[,i])#提取每个变量的最小值
    for (j in 1:nrow)
    {
      new[j,i]=(x[j,i]-min)/(max-min)#计算归一化后的新数据集
    }
  }
  new
}
wine=na.omit(data)
wine<-subset(wine, select = c("id","TRIP10","NGFR","SLC48A1","SRMS","riskScore","risk"))
row.names(wine)<-wine$id
wine<-wine[,-1]

set.seed(71)
samp=sample(1:518,518)
#wine=wine[,-1]
#wine=subset(wine,select = c("futime","TRIP10","NGFR","SLC48A1","SRMS","riskScore","V11"))
wine[samp,1:5]=scale01(wine[samp,])#对样本进行预处理
r=1/max(abs(wine[samp,1:5]))#确定参数rang的变化范围
set.seed(101)
model1=nnet(risk~.,data=wine,subest=samp,size=4,rang=r,decay=5e-4,maxit=200)#建立神经网络模型

#提取wine数据集中除quality列以外的数据作为自变量
x=subset(wine,select=-risk)
y=wine[,6]#提取quality列数据作为响应变量
y=class.ind(y)#预处理将其变为类指标矩阵
set.seed(101)
model2=nnet(x,y,decay=5e-4,maxit=200,size=4,rang=r)#建立神经网络模型

summary(model1)
summary(model2)


#预测判别
#第一种建模方式建立的模型
x=wine[,1:5]
pred=predict(model1,x,type="class")#根据模型对x数据进行预测
set.seed(110)
pred[sample(1:518,8)]

#第二种建模方式建立的模型
xt=wine[,1:5]
pred=predict(model2,xt)#根据模型对xt数据进行预测
dim(pred)
pred[sample(1:518,4),]
name=c("high","low")
prednew=max.col(pred)#确定每行中最大值所在的列
prednewn=name[prednew]#根据预测结果将其变为相对应的类别名称
set.seed(201)
prednewn[sample(1:518,8)]
true=max.col(y)#确定真实值的每行中最大值所在的列
table(true,prednewn)#模型预测精度展示

#优化模型
#size=i控制隐藏层节点
set.seed(444)
nrow.wine=dim(wine)[1]
samp=sample(1:nrow.wine,nrow.wine*0.7)#抽取70%样本
wine[samp,1:5]=scale01(wine[samp,])#对数据样本进行预处理
wine[-samp,1:5]=scale01(wine[-samp,])#对测试集进行预处理
r=1/max(abs(wine[samp,1:5]))#确定rang的变化范围
n=length(samp)
err1=0
err2=0
for(i in 1:17)
{
  set.seed(111)
  model=nnet(risk~.,data=wine,maxit=400,rang=r,size=i,subset=samp,decay=5e-4)
  err1[i]=sum(predict(model,wine[samp,1:5],type='class')!=wine[samp,6])/n
  err2[i]=sum(predict(model,wine[-samp,1:5],type='class')!=wine[-samp,6])/(nrow.wine-n)
}  #运行时间较长，

plot(1:17,err1,'l',col=1,lty=1,ylab="模型误判率",xlab="隐藏层节点个",ylim=c(min(min(err1),min(err2)),max(max(err1),max(err2))))
lines(1:17,err2,col=1,lty=3)
points(1:17,err1,col=1,pch="+")
points(1:17,err2,col=1,pch="o")
legend(1,0.53,"测试集误判率",bty="n",cex=1.5)
legend(1,0.35,"训练集误判率",bty="n",cex=1.5)


# maxit:控制的是 模型的最大迭代次数
err11=0
err12=0
for(i in 1:500)
{
  set.seed(111)   
  model=nnet(risk~.,data=wine,maxit=i,rang=r,size=8,subset=samp)
  err11[i]=sum(predict(model,wine[samp,1:5],type='class')!=wine[samp,6])/n
  err12[i]=sum(predict(model,wine[-samp,1:5],type='class')!=wine[-samp,6])/(nrow.wine-n)
}

plot(1:length(err11),err11,'l',ylab="模型误判率",xlab="训练周期",col=1,ylim=c(min(min(err11),min(err12)),max(max(err11),max(err12))))
lines(1:length(err11),err12,col=1,lty=3)
legend(250,0.50,"测试集误判率",bty="n",cex=1.2)
legend(250,0.45,"训练集误判率",bty="n",cex=1.2)


# 取迭代次数为20，隐藏节点为8
set.seed(111)
model=nnet(risk~.,data=wine,maxit=20,rang=r,size=8,subset=samp)
x=wine[-samp,1:5]
pred=predict(model,x,type="class")
table(wine[-samp,6],pred)

mydata = data.frame(
  prov = c("nobenefit", "benefit"),
  sensitive = c(22/78,78/78),
  resistant = c(56/78, 0/78)
)
# 获取数据结构
str(mydata)
library("ggplot2")
library("tidyr")
boxplot_data <- gather(mydata, E1, E2, -prov) 
pdf("ANN.boxplot.pdf")
ggplot(boxplot_data) +
  geom_bar(aes(x = prov, y = E2, fill = E1),
           stat = "identity") +
  theme_bw() 
dev.off()

testing<-wine[-samp,]
table(testing$risk)
#high  low 
#78   78

testing$prediction<-pred
testing$id<-row.names(testing)
testing_surv<-merge(testing,data,by="id")

library("survival")
library("survminer")
risk<-subset(testing_surv,select = c("futime", "fustat", "prediction"))
rt=risk
diff=survdiff(Surv(futime, fustat) ~prediction,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(futime, fustat) ~ prediction, data = rt)

#绘制生存曲线
pdf(file="ANN.survival.pdf",onefile = FALSE,
    width = 5.5,             #图片的宽度
    height =5)             #图片的高度
ggsurvplot(fit, 
           data=rt,
           conf.int=TRUE,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=TRUE,
           legend.labs=c("No benefit", "Benefit"),
           legend.title="Risk",
           xlab="Time(years)",
           break.time.by = 1,
           risk.table.title="",
           palette=c("red", "blue"),
           risk.table.height=.25)
dev.off()

#############################################################################################
exterdata<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/05_C-index/GSE143985_GPL570/expData_validation_1.csv",header = T, row.names = 1)
riskgene<-c("TRIP10","NGFR","SLC48A1","SRMS")
extdata<-exterdata[riskgene,]
extdata_t<-t(extdata)
extdata_t<-as.data.frame(extdata_t)
extdata_t$id<-row.names(extdata_t)
ext_cli<-read.table("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/06_Time-dependent_ROC/risk_GSE143985.txt",header=T)
ext_merg<-merge(extdata_t,ext_cli,by="id")
ext_merg$risk<-factor(ext_merg$risk)
ext_merg<-subset(ext_merg,select = c("id","TRIP10","NGFR","SLC48A1.x","SRMS","dfs_time","dfs_event","riskScore","risk"))
colnames(ext_merg)<-c("id","TRIP10","NGFR","SLC48A1","SRMS","dfs_time","dfs_event","riskScore","risk")
ext_pre<-subset(ext_merg, select = c("id","TRIP10","NGFR","SLC48A1","SRMS","riskScore","risk"))
row.names(ext_pre)<-ext_pre$id
ext_pre<-ext_pre[,-1]

set.seed(71)
samp=sample(1:91,91)
ext_pre[samp,1:5]=scale01(ext_pre[samp,])#对样本进行预处理
r=1/max(abs(ext_pre[samp,1:5]))#确定参数rang的变化范围
set.seed(101)
model1=nnet(risk~.,data=ext_pre,subest=samp,size=4,rang=r,decay=5e-4,maxit=200)#建立神经网络模型


x=subset(ext_pre,select=-risk)
y=ext_pre[,6]#提取quality列数据作为响应变量
y=class.ind(y)#预处理将其变为类指标矩阵
set.seed(101)
model2=nnet(x,y,decay=5e-4,maxit=200,size=4,rang=r)#建立神经网络模型

summary(model1)
summary(model2)


#预测判别
#第一种建模方式建立的模型
x=ext_pre[,1:5]
pred=predict(model1,x,type="class")#根据模型对x数据进行预测
set.seed(110)
pred[sample(1:91,8)]


#第二种建模方式建立的模型
xt=ext_pre[,1:5]
pred=predict(model2,xt)#根据模型对xt数据进行预测
dim(pred)
pred[sample(1:91,4),]
name=c("high","low")
prednew=max.col(pred)#确定每行中最大值所在的列
prednewn=name[prednew]#根据预测结果将其变为相对应的类别名称
set.seed(201)
prednewn[sample(1:91,8)]
true=max.col(y)#确定真实值的每行中最大值所在的列
table(true,prednewn)#模型预测精度展示


#优化模型
#size=i控制隐藏层节点
set.seed(444)
nrow.ext_pre=dim(ext_pre)[1]
samp=sample(1:nrow.ext_pre,nrow.ext_pre*0.7)#抽取70%样本
ext_pre[samp,1:5]=scale01(ext_pre[samp,])#对数据样本进行预处理
ext_pre[-samp,1:5]=scale01(ext_pre[-samp,])#对测试集进行预处理
r=1/max(abs(ext_pre[samp,1:5]))#确定rang的变化范围
n=length(samp)
err1=0
err2=0
for(i in 1:17)
{
  set.seed(111)
  model=nnet(risk~.,data=ext_pre,maxit=400,rang=r,size=i,subset=samp,decay=5e-4)
  err1[i]=sum(predict(model,ext_pre[samp,1:5],type='class')!=ext_pre[samp,6])/n
  err2[i]=sum(predict(model,ext_pre[-samp,1:5],type='class')!=ext_pre[-samp,6])/(nrow.ext_pre-n)
}  #运行时间较长，

plot(1:17,err1,'l',col=1,lty=1,ylab="模型误判率",xlab="隐藏层节点个",ylim=c(min(min(err1),min(err2)),max(max(err1),max(err2))))
lines(1:17,err2,col=1,lty=3)
points(1:17,err1,col=1,pch="+")
points(1:17,err2,col=1,pch="o")
legend(1,0.53,"测试集误判率",bty="n",cex=1.5)
legend(1,0.35,"训练集误判率",bty="n",cex=1.5)


# maxit:控制的是 模型的最大迭代次数
err11=0
err12=0
for(i in 1:500)
{
  set.seed(111)   
  model=nnet(risk~.,data=ext_pre,maxit=i,rang=r,size=8,subset=samp)
  err11[i]=sum(predict(model,ext_pre[samp,1:5],type='class')!=ext_pre[samp,6])/n
  err12[i]=sum(predict(model,ext_pre[-samp,1:5],type='class')!=ext_pre[-samp,6])/(nrow.ext_pre-n)
}

plot(1:length(err11),err11,'l',ylab="模型误判率",xlab="训练周期",col=1,ylim=c(min(min(err11),min(err12)),max(max(err11),max(err12))))
lines(1:length(err11),err12,col=1,lty=3)
legend(250,0.50,"测试集误判率",bty="n",cex=1.2)
legend(250,0.45,"训练集误判率",bty="n",cex=1.2)


# 取迭代次数为20，隐藏节点为1
set.seed(111)
model=nnet(risk~.,data=ext_pre,maxit=20,rang=r,size=8,subset=samp)
x=ext_pre[,1:5]
pred=predict(model,x,type="class")
table(ext_pre[,6],pred)

mydata = data.frame(
  prov = c("nobenefit", "benefit"),
  sensitive = c(36/45,46/46),
  resistant = c(9/45,0/46)
)
# 获取数据结构
str(mydata)
library("ggplot2")
library("tidyr")
boxplot_data <- gather(mydata, E1, E2, -prov) 
pdf("ANN.boxplot_ext.pdf")
ggplot(boxplot_data) +
  geom_bar(aes(x = prov, y = E2, fill = E1),
           stat = "identity") +
  theme_bw() 
dev.off()

ext_pre$prediction<-pred
ext_pre$id<-row.names(ext_pre)
ext_pre_surv<-merge(ext_pre,ext_merg,by="id")
 
risk<-subset(ext_pre_surv,select = c("dfs_time", "dfs_event", "prediction"))
rt=risk
diff=survdiff(Surv(dfs_time, dfs_event) ~prediction,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)

fit <- survfit(Surv(dfs_time, dfs_event) ~ prediction, data = rt)

#绘制生存曲线
pdf(file="ANN.survival_ext.pdf",onefile = FALSE,
    width = 5.5,             #图片的宽度
    height =5)             #图片的高度
ggsurvplot(fit, 
           data=rt,
           conf.int=TRUE,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=TRUE,
           legend.labs=c("No benefit", "Benefit"),
           legend.title="Risk",
           xlab="Time(years)",
           break.time.by = 1,
           risk.table.title="",
           palette=c("red", "blue"),
           risk.table.height=.25)
dev.off()



