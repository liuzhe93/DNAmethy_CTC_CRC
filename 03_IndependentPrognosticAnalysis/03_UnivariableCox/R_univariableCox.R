setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/03_IndependentPrognosticAnalysis/03_UnivariableCox")
rm(list=ls())

library(survival)
library(forestplot)

clrs <- fpColors(box="green",line="darkblue", summary="royalblue")             #定义森林图颜色
rt=read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/03_IndependentPrognosticAnalysis/02_Time-dependent_ROC/indepInput.csv",header=T)
head(rt)
for(i in 1:(dim(rt)[1])){
  rt$sample_name_new[i]<-paste(rt$sample_name[i],i,sep="_")
}
rownames(rt)<-rt$sample_name_new
rt<-subset(rt, select = c("futime", "Status", "Age", "Gender", "Stage", "riskScore"))
colnames(rt)<-c("futime", "fustat", "Age", "Gender", "Stage", "riskScore")
head(rt)

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(outTab,file="uniCox.xls",sep="\t",row.names=F,quote=F)

#绘制森林图
rt=read.table("uniCox.xls",header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(rt)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )          #定义图片文字
pdf(file="forest.pdf",onefile = FALSE,
    width = 6,             #图片的宽度
    height = 4,            #图片的高度
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()
