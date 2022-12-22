setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/05_C-index")
rm(list=ls())


###############################################TCGA C-index calculation##########################################
com_auc<-read.csv("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/03_Gaussian_mixture_model-based_hierarchical_clustering_method/cluster_results.csv",
                  header=T)
com_auc_sorted=com_auc[order(com_auc$auc,decreasing = T),]
genes<-com_auc_sorted$gene[1]
gene_list<-strsplit(genes,split="\\ \\+\\ ")[[1]]
uniSigExp<-read.table("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/01_univariate_Cox_regression_analysis/uniSigExp.txt",
                      sep="\t",header=T)
uniSigExp_filtered<-uniSigExp[,c("id","futime","fustat",gene_list)]
dim(uniSigExp_filtered)
#[1] 521   7
head(uniSigExp_filtered)
library("survival")
fit <- coxph(Surv(futime, fustat) ~ TRIP10 + NGFR + SLC48A1 + SRMS, data = uniSigExp_filtered)
sum.surv <- summary(fit)
c_index <- sum.surv$concordance
c_index
#         C      se(C) 
#0.66920668 0.02680977
############################validation dataset 1: GSE143985_GPL570#######################################################
setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/05_C-index/GSE143985_GPL570")
clinical<-read.csv("clinical.csv", header = F)
clinical_infor<-t(clinical)
colnames(clinical_infor)<-clinical_infor[1,]
clinical_infor<-clinical_infor[-1,]
clinical_infor<-as.data.frame(clinical_infor)
clinical_infor$dfs_time<-as.numeric(clinical_infor$dfs_time)/365
write.csv(clinical_infor, "clinical_infor.csv", quote=F, row.names = F)

library("affy")
data <- ReadAffy(celfile.path="GSE143985_RAW")
est <- rma(data)
exprSet <- exprs(est)
dim(exprSet)
#[1] 54675    91
#phenoData(est)
#sampleNames(est)
#varMetadata(est) 
featureNames(est)[1:10]
library("dplyr")
featureNames(est) %>% length() 
#[1] 54675
featureNames(est) %>% unique() %>% length() 
#[1] 54675
# 根据不同的平台号选择相应的R包
# http://www.bio-info-trainee.com/1399.html
library("hgu133plus2.db")
ls("package:hgu133plus2.db")
ids=toTable(hgu133plus2SYMBOL)
head(ids)
length(unique(ids$symbol))
#[1] 20857
tail(sort(table(ids$symbol)))
#每个基因对应多少探针
#  ZBTB20 ARHGEF12     CD44    CFLAR    DNAH1      HFE 
#      12       13       13       13       13       15
table(sort(table(ids$symbol)))
#   1    2    3    4    5    6    7    8    9   10   11   12   13   15 
#9815 5310 2879 1434  755  360  171   66   34   15    9    4    4    1
#9815个基因设计了1个探针；5310个基因设计了2个探针；2879个基因设计了3个探针.....也就是说大部分的基因只设计了1个探针
table(rownames(exprSet) %in% ids$probe_id)
#FALSE  TRUE 
#11537 43138
#发现有11537个探针没有对应的基因名；43138个探针有对应的基因名。
#现在我们对探针进行过滤，把没有对应基因名的探针过滤掉：
exprSet_filter = exprSet[rownames(exprSet) %in% ids$probe_id,]
dim(exprSet)
#[1] 54675    91
dim(exprSet_filter)
#[1] 43138    91
ids=ids[match(rownames(exprSet_filter),ids$probe_id),]
#然后，我们使用match函数把ids里的探针顺序改一下，使ids里探针顺序和我们表达矩阵的顺序完全一样。
head(ids)
head(exprSet_filter)
colnames(exprSet_filter)<-substr(colnames(exprSet_filter),1,10)
exprSet_filter<-as.data.frame(exprSet_filter)
tmp = by(exprSet_filter,
         ids$symbol,
         function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character(tmp)
dim(exprSet_filter)
exprSet_filter = exprSet_filter[rownames(exprSet_filter) %in% probes,] # 过滤有多个探针的基因
dim(exprSet_filter)
#[1] 20857    91
head(tmp)
head(probes)
rownames(exprSet_filter)=ids[match(rownames(exprSet_filter),ids$probe_id),2]
write.csv(exprSet_filter, file="expData_validation_1.csv",quote=F)
exp_selected<-exprSet_filter[gene_list,]
exp_selected_t<-t(exp_selected)
exp_selected_t<-as.data.frame(exp_selected_t)
exp_selected_t$sample_name<-row.names(exp_selected_t)
mydata_vali_1<-merge(clinical_infor,exp_selected_t,by="sample_name")
mydata_vali_1$dfs_event<-as.integer(mydata_vali_1$dfs_event)
library("survival")
fit <- coxph(Surv(dfs_time, dfs_event) ~ TRIP10 + NGFR + SLC48A1 + SRMS, data = mydata_vali_1)
sum.surv <- summary(fit)
c_index <- sum.surv$concordance
c_index
#         C      se(C) 
#0.77284372 0.05681457
############################validation dataset 2: GSE63624_GPL5175#######################################################
setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/05_C-index/GSE63624_GPL5175")
clinical<-read.csv("clinical.csv", header = F)
clinical_infor<-t(clinical)
colnames(clinical_infor)<-clinical_infor[1,]
clinical_infor<-clinical_infor[-1,]
clinical_infor<-as.data.frame(clinical_infor)
clinical_infor$rfs_time<-as.numeric(clinical_infor$rfs_time)/365
write.csv(clinical_infor, "clinical_infor.csv", quote=F, row.names = F)

library("oligo")
data.dir<-"GSE63624_RAW"
celfiles<-list.files(data.dir,"\\.gz$")
data.raw<-read.celfiles(filenames = file.path(data.dir,celfiles))
sampleNames(data.raw)<-substr(sampleNames(data.raw),1,10)
#pData(data.raw)
data.eset<-rma(data.raw)
data.exprs<-exprs(data.eset)
dim(data.exprs)
#[1] 22011    52
str(data.exprs)
head(data.exprs)
#P/A过滤：目的是去除“不表达”的基因/探针数据，使用paCalls函数，选取p值小于0.05的探针：
xpa <- paCalls(data.raw)
head(xpa)
AP <- apply(xpa, 1, function(x) any(x < 0.05))
xids <- as.numeric(names(AP[AP]))
head(xids)
pinfo<-getProbeInfo(data.raw)
head(pinfo)
#两类id转换后进行筛选：
fids <- pinfo[pinfo$fid %in% xids, 2]
head(fids)
nrow(data.exprs)
#[1] 22011
data.exprs <- data.exprs[rownames(data.exprs) %in% fids, ]
nrow(data.exprs)
#[1] 21838
write.csv(data.exprs, file="expData_validation_2.csv",quote=F)
#data.exprs<-read.csv("expData_validation_2.csv",header=T,row.names = 1)

# 探针到基因的转换
gpl_infor<-read.table("GPL5175-3188.txt",header=T,comment.char = "#",sep="\t",fill=T)
ids = gpl_infor[,c("ID","gene_assignment")]
gene_list
#[1] "TRIP10"  "NGFR"    "SLC48A1" "SRMS" 
#"TRIP10"
#3818515

#"NGFR"
#3725685

#"SLC48A1"
#3413212

#"SRMS" 
#3913945 
TRIP10_exp<-data.exprs["3818515",]
NGFR_exp<-data.exprs["3725685",]
SLC48A1_exp<-data.exprs["3413212",]
SRMS_exp<-data.exprs["3913945",]
exp_selected<-rbind(TRIP10_exp,NGFR_exp,SLC48A1_exp,SRMS_exp)
exp_selected_t<-t(exp_selected)
exp_selected_t<-as.data.frame(exp_selected_t)
colnames(exp_selected_t)<-c("TRIP10","NGFR","SLC48A1","SRMS")
exp_selected_t$sample_name<-row.names(exp_selected_t)
mydata_vali_2<-merge(clinical_infor,exp_selected_t,by="sample_name")
mydata_vali_2$rfs_event<-as.integer(mydata_vali_2$rfs_event)
library("survival")
fit <- coxph(Surv(rfs_time, rfs_event) ~ TRIP10 + NGFR + SLC48A1 + SRMS, data = mydata_vali_2)
sum.surv <- summary(fit)
c_index <- sum.surv$concordance
c_index
#        C      se(C) 
#0.58392435 0.05958164

########################################################################################################
setwd("D:/MyProjects/CTC/CTC_Prognosis_COAD/02_ConstructionRiskModel/05_C-index")
table = data.frame(
  tcga = c("TCGA", 0.66920668),
  gse1 = c("GSE143985", 0.77284372)
#  gse2 = c("GSE63624",0.58392435)
)
print(table) # 查看 table 数据
pdf("C-index.pdf")
barplot(height = c(0.66920668, 0.77284372),  # 绘图数据（矩阵）
        names.arg = c('TCGA', 'GSE143985'),  # 柱子名称
        col = 'steelblue',  # 填充颜色
        border = '#ffffff',   # 轮廓颜色
        xlab = 'dataset',  # X轴名称
        ylab = 'C-index',  # Y轴名称
        main = 'C-index',  # 主标题
        horiz = FALSE,  # 是否为水平放置
        ylim = c(0, 1), # Y轴取值范围
)
dev.off()

