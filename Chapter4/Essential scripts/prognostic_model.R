setwd('/Users/zhuofanmou/Documents/PhD Exeter/ClariomD/alternative_splicing/R_workspace_v2/survival_analysis')
getwd()

# read in prepared Clinical info
clinic_info_AS = read.table("clinical_OS.txt",header = T,sep="\t",row.names = 1,check.names = F)
View((clinic_info_AS))

# Preparation 
clinic_info_AS <- clinic_info_AS[ clinic_info_AS$shortLetterCode !='NT', ] # 499 with metastatic tumor included
clinic_info_AS <- clinic_info_AS[ clinic_info_AS$shortLetterCode !='TM', ] # 498 with metastatic tumor excluded
dim(clinic_info_AS)

# NOTE we have 3 patients with replicated tumour samples (A&B), so we need to keep one (A), so 498=>495
clinic_info_AS <- clinic_info_AS[!duplicated(clinic_info_AS$patient),]

#clinic_info_AS$sampleID <- rownames(clinic_info_AS) # we need this so later we can merge with asMatrix from 

rownames(clinic_info_AS) <- clinic_info_AS$patient


clinic_info_AS[clinic_info_AS$vital_status == 'Alive',]$vital_status <- 0 
clinic_info_AS[clinic_info_AS$vital_status == 'Dead',]$vital_status <- 1
table(clinic_info_AS$vital_status)
#clinic_info_AS <- clinic_info_AS[ clinic_info_AS$vital_status !='Not Reported', ] # 493 with 2 patient without OS vital status reported 'Not Reported'; BUT we dont need to remove this if we do for the progression-free survivsl
dim(clinic_info_AS)
str(clinic_info_AS)

# add PFS 
load_pfs <- load("TCGA_PAN_PFS_GSig_Input_edgeR.rda")
load_pfs

pfs <- rt_final_PFS_edgeR

View(pfs)
pfs <- pfs[ , -c(3:12)]
pfs$name <- rownames(pfs)
pfs$name = gsub(pattern = "-01.*", replacement = "", x = pfs$name)
rownames(pfs) <- pfs$name
pfs <- pfs[ , -3]

#clinic_info_AS$vital_status <- as.numeric(clinic_info_AS$vital_status)

colnames(pfs)[1] <- ('futime')
colnames(pfs)[2] <- ('fustat')



##
sameSample_AS=intersect(row.names(clinic_info_AS),row.names(pfs)) 
length(sameSample_AS)

##

##
df1=clinic_info_AS[sameSample_AS,] 
dim(df1)
head(rownames(df1))
df2 <- pfs[sameSample_AS,] 
dim(df2)
head(rownames(df2))
##


##Check if row names matched
table(row.names(df1) == row.names(df2))
##

## Combine the data (including PSA) for each of the cohorts
df_clinicPFS_AS <- cbind(df1,df2) 
View(df_clinicPFS_AS)
dim(df_clinicPFS_AS)
##

df_clinicPFS_AS <- df_clinicPFS_AS[ , -c(1:5)]
df_clinicPFS_AS <- df_clinicPFS_AS[ , c(8,9,1:7)]

write.table(df_clinicPFS_AS,
            file = "clinical_final_PFS_forAS.txt",
            sep = "\t",
            quote = F,
            row.names = T) 


# PFS, BCRFS and OS will be combined in 'differentDataset_clinicalCompare.R'


#==== Import TCGA SpliceSeq pre-processed full-AS-PSI, and common-AS-PSI expression matrix between TCGAss and PRJEB2449

#8440 DEAS events identified from 'DEAS.R'
DEAS=read.delim("DEAS.txt",sep="\t",header=T,check.names=F,row.names=1) 
View(DEAS)
dim(DEAS)


# all 29415 events PSI and 541 sampels 
all=read.table("ASevents_PSIdata.txt",sep="\t",header=T,check.names=F,row.names=1)    
all = subset(all, rownames(all) %in% rownames(DEAS)) # subsetting all events to the 8440 DEAS events
dim(all)

# matched events of DEAS between TCGAss and PRJEB2449
matched_Events <- read.delim("CommonEvents_in_All.txt")
View(matched_Events)
rownames(matched_Events) <- matched_Events$ID
all = subset(all, rownames(all) %in% rownames(matched_Events)) # subsetting DEAS events to the 141 matched events
dim(all)
all = rbind(id = colnames(all), all)
View(all)

##################################################################
########################### survival data  #######################
##################################################################
#=== 141 PSI matrix
View(all)
PSIexp <- all

PSIexp <- t(PSIexp)
dim(PSIexp)
View(PSIexp)
PSI_TCGA <- save(PSIexp, file = "PSI_TCGA.rda")

#=== import all clinical data (pre-prepared)
survival_data<- read.delim("PFSpaper1vsTCGAbiolinkslimour_allClinicopathoVars.txt")
View(survival_data)

#=== Progreaaion-free survival (PFS) and combine PFS with 141 PSI matrix
PFS_df <- survival_data[ , c(1,2)]
View(PFS_df)

## Combine
sameSample_PFS=intersect(row.names(PFS_df),row.names(PSIexp)) 
length(sameSample_PFS)

df1_pfs=PFS_df[sameSample_PFS,] 
dim(df1_pfs)
head(rownames(df1_pfs))
df2_pfs <- PSIexp[sameSample_PFS,] 
dim(df2_pfs)
head(rownames(df2_pfs))

##Check if row names matched
table(row.names(df1_pfs) == row.names(df2_pfs))

## Combine the data (including PSA) for each of the cohorts
asTime_PFS <- cbind(df1_pfs,df2_pfs) 
View(asTime_PFS)
asTime_PFS <- asTime_PFS[ , c(3,1,2,4:ncol(asTime_PFS))]
colnames(asTime_PFS)[1] <- ('id')   # the complete PFS combined with PSI values data matrix 
View(head(asTime_PFS))
dim(asTime_PFS) # 413*144
table(asTime_PFS$fustat) # 84

# save PSI matrix with PFS information (the complete set: 413)
write.table(asTime_PFS,
            file = "asTime_PFS.txt",
            sep = "\t",
            quote = F,
            row.names = F) 


#=== Biochemical recurrence free survival (BCRFS) and combine BCRFS with 141 PSI matrix
BCRFS_df <- survival_data[ , c(24,23)]
View(BCRFS_df)
BCRFS_df$bcrf_time.1 <- BCRFS_df$bcrf_time.1/365
colnames(BCRFS_df)[1] <- ('futime')
colnames(BCRFS_df)[2] <- ('fustat')
table(BCRFS_df$fustat)

## Combine
sameSample_BCRFS=intersect(row.names(BCRFS_df),row.names(PSIexp)) 
length(sameSample_BCRFS)

df1_bcrfs=BCRFS_df[sameSample_BCRFS,] 
dim(df1_bcrfs)
head(rownames(df1_bcrfs))
df2_bcrfs <- PSIexp[sameSample_BCRFS,] 
dim(df2_bcrfs)
head(rownames(df2_bcrfs))

##Check if row names matched
table(row.names(df1_bcrfs) == row.names(df2_bcrfs))

## Combine the data (including PSA) for each of the cohorts
asTime_BCRFS <- cbind(df1_bcrfs,df2_bcrfs) 
View(asTime_BCRFS)
asTime_BCRFS <- asTime_BCRFS[ , c(3,1,2,4:ncol(asTime_BCRFS))]
colnames(asTime_BCRFS)[1] <- ('id')   # the complete BCRFS combined with PSI values data matrix 
View(head(asTime_BCRFS))
dim(asTime_BCRFS) # 413*144
table(asTime_BCRFS$fustat) # 67

# save PSI matrix with BCRFS information (the complete set: 413)
write.table(asTime_BCRFS,
            file = "asTime_BCRFS.txt",
            sep = "\t",
            quote = F,
            row.names = F) 


#=== Overall survival (OS) and combine OS with 141 PSI matrix
OS_df <- survival_data[ , c(20,19)]
View(OS_df)
OS_df$os_time.1 <- OS_df$os_time.1/365
colnames(OS_df)[1] <- ('futime')
colnames(OS_df)[2] <- ('fustat')
table(OS_df$fustat)


## Combine
sameSample_OS=intersect(row.names(OS_df),row.names(PSIexp)) 
length(sameSample_OS)

df1_os=OS_df[sameSample_OS,] 
dim(df1_os)
head(rownames(df1_os))
df2_os <- PSIexp[sameSample_OS,] 
dim(df2_os)
head(rownames(df2_os))

##Check if row names matched
table(row.names(df1_os) == row.names(df2_os))

## Combine the data (including PSA) for each of the cohorts
asTime_OS <- cbind(df1_os,df2_os) 
View(asTime_OS)
asTime_OS <- asTime_OS[ , c(3,1,2,4:ncol(asTime_OS))]
colnames(asTime_OS)[1] <- ('id')   # the complete BCRFS combined with PSI values data matrix 
View(head(asTime_OS))
dim(asTime_OS) # 413*144
table(asTime_OS$fustat) # 9

# save PSI matrix with OS information (the complete set: 413)
write.table(asTime_OS,
            file = "asTime_OS.txt",
            sep = "\t",
            quote = F,
            row.names = F) 


############################################### This script FINISHES HERE!!! ##################################################


#################################################################################################################
##################### UniCox regression (immAS10.uniCox.R) ######################################################
################################################################################################################


# library
library(survival)
library(UpSetR)
inputFile="asTime.txt"
pFilter=0.05

# read in data
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
View(rt)
#Univariate Cox regression
outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  cox=coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
  coxSummary=summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     z=coxSummary$coefficients[,"z"],
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}

#output from Univariate Cox regression analysis
outTab=outTab[is.na(outTab$pvalue)==FALSE,]
outTab=outTab[order(as.numeric(as.vector(outTab$pvalue))),]
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)
#output for significant events from Uni Cox regression analysis
sigTab=outTab[as.numeric(as.vector(outTab$pvalue))<pFilter,]
write.table(sigTab,file="uniCox.Sig.txt",sep="\t",row.names=F,quote=F)
# Uni-Cox significant AS events with their PSI
sigGenes=c("futime", "fustat", as.vector(sigTab[,1]))
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp), uniSigExp)
write.table(uniSigExp, file="uniSigExp.txt", sep="\t", row.names=F, quote=F)


# UpeSet for Uni-Cox significant events
gene=sapply(strsplit(sigGenes,"\\|"),"[",1)
asType=sapply(strsplit(sigGenes,"\\|"),"[",3)
upsetList=list(AA=unique(gene[asType=="AA"]),
               AD=unique(gene[asType=="AD"]),
               AP=unique(gene[asType=="AP"]),
               AT=unique(gene[asType=="AT"]),
               ES=unique(gene[asType=="ES"]),
               ME=unique(gene[asType=="ME"]),
               RI=unique(gene[asType=="RI"]) )
upsetData=fromList(upsetList)

#输出图形
pdf(file="uniCoxUpset.pdf", width=8, height=5, onefile=FALSE)
upset(upsetData,
      nsets = 7,                    #展示可变剪切类型个数
      order.by = "freq",            #按照基因数目排序
      show.numbers = "yes",         #柱状图上方是否显示数值
      number.angles = 20,           #字体角度
      point.size = 1.5,             #点的大小
      matrix.color="red",           #交集点颜色
      line.size = 0.8,              #线条粗线
      mainbar.y.label = "Gene Intersections",
      sets.x.label = "Set Size")
dev.off()



################################################################################################################
########################## bubble plot (immAS11.bubble.R) ######################################################
################################################################################################################


# library(ggplot2)           #引用包
# inputFile="uniCox.txt"     #单因素结果文件
# pFilter=0.05               #单因素显著性过滤条件
# #setwd("C:\\biowolf\\immAS\\11.bubble")            #设置工作目录
# rt=read.table(inputFile, header=T, sep="\t", check.names=F)     #读取输入文件
# 
# #定义显著性
# Significant=ifelse(rt$pvalue<pFilter, "Prognosis AS", "No significant")
# #绘制火山图
# p=ggplot(rt, aes(z, -log10(pvalue)))+
#   geom_point(aes(col=Significant))+
#   scale_color_manual(values=c("#2f5688","#CC0000"))+
#   labs(title=" ")+ xlab("z-score")+
#   geom_vline(xintercept=0, linetype="dotted")+
#   theme(plot.title = element_text(size=16, hjust=0.5, face="bold"))
# p=p+theme_bw()
# #保存为图片
# pdf("vol.pdf", width=6, height=5)
# print(p)
# dev.off()
# 
# 
# #准备气泡图输入数据
# rt=rt[order(as.numeric(as.vector(rt$pvalue))),]
# rt=rt[rt$pvalue<pFilter,]
# row.names(rt)=rt[,1]
# #rt[,1]=gsub("\\|", "\\-", rt[,1])
# asTypes=gsub("(.*)\\|(.*)\\|(.*?)", "\\3", row.names(rt)) 
# 
# #绘制气泡图
# for(asType in levels(factor(asTypes)) ){
#   genes=rownames(rt)
#   gene=grep(paste0("\\|",asType), genes, value=T)
#   geneLength=ifelse(length(gene)>20, 20, length(gene))
#   data=rt[gene[1:geneLength],]
#   
#   data=data[order(as.numeric(as.vector(data$pvalue)),decreasing=T),]
#   data$id = factor(data$id,levels=as.character(data[,1]))
#   
#   #绘制气泡图
#   p=ggplot(data,aes(z,id))		
#   pbubble = p + geom_point(aes(color=pvalue, size=-1*log10(pvalue)) )
#   pr = pbubble + 
#     scale_colour_gradient(low="red", high="skyblue") + 
#     labs(color="pvalue", size="-log10(pvalue)", x="z-score", y="")+
#     guides(color = guide_colourbar(order = 1), size = guide_legend(order = 2))+
#     theme_bw()
#   #输出图形
#   ggsave(paste0("bubble.", asType, ".pdf"),width=5.5,height=5)
# }



################################################################################################################
################################# model (immAS12.model.R) ######################################################
################################################################################################################

#引用包
library("glmnet")
library("survival")
asType=""                     #(AA AD AP AT ES RI ME)中的一种,如果7种一起做，asType设置为空即可
inputFile= "asTime.txt"     #"uniSigExp.txt"     #输入文件
#setwd("C:\\biowolf\\immAS\\12.model")      #设置工作目录

#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
View(rt)
# rt$futime[rt$futime<=0]=1
# rt$futime=rt$futime/365
genes=colnames(rt)
gene=grep(paste0("\\|",asType), genes, value=T)
geneLength=ifelse(length(gene)>20, 20, length(gene))
rt=rt[,c("futime", "fustat", gene[1:geneLength])]

rt <- rt[ , c("futime","fustat","TTC31|54078|RI",
        "PSMD8|49638|AT",
        "PDCD2|78501|RI",
        "CDC42BPA|10042|ES",
        "TGIF1|44505|AT",
        "SLC25A19|43430|AP",
        "TTLL5|28520|AT",
        "TMUB2|41787|AT")]

# set a seed so lasso will be the same and reviewer can retrieve
set.seed(123)
#lasso回归
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit <- glmnet(x, y, family = "cox", maxit = 1000)

#lasso图形
pdf("lasso.lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

#交叉验证图形
cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("lasso.cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

#输出lasso显著基因的表达文件
coef <- coef(fit, s = cvfit$lambda.min) # use cvfit will be the same, why???
coef
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene=c("futime","fustat",lassoGene)
lassoSigExp=rt[,lassoGene]
lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
View(lassoSigExp)
write.table(lassoSigExp,file="lasso.sigExp.txt",sep="\t",row.names=F,quote=F)

#COX模型构建
rt=read.table("lasso.sigExp.txt",header=T,sep="\t",check.names=F,row.names=1)
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)#coxph(Surv(futime, fustat) ~ `PDCD2|78501|RI`, data = rt)#coxph(Surv(futime, fustat) ~ ., data = rt)


# # stepwise regression: https://www.jianshu.com/p/33425451f1f2?utm_campaign=hugo&utm_medium=reader_share&utm_content=note
# 
# colnames(rt) <- c("futime","fustat","TTC31.54078.RI",
#               "PSMD8.49638.AT",
#               "PDCD2.78501.RI",
#               "CDC42BPA.10042.ES",
#               "TGIF1.44505.AT")
# 
# list <- colnames(rt)[c(3:ncol(rt))]
# 
# library(My.stepwise)
# My.stepwise.coxph(Time = "futime", Status = "fustat", variable.list = list, data = rt)
# # # ================================================================================================== 
# # *** Stepwise Final Model (in.lr.test: sle = 0.15; out.lr.test: sls = 0.15; variable selection restrict in vif = 999): 
# #   Call:
# #   coxph(formula = Surv(futime, fustat) ~ PDCD2.78501.RI + CDC42BPA.10042.ES, 
# #         data = data, method = "efron")

multiCox=step(multiCox,direction = "both")

# ## Use LASSO-COX regression with FIX coefficients!!! https://stats.stackexchange.com/questions/572211/questions-related-to-survival-analysis-and-lasso-cox-regression
# coefs = coef(cvfit, s = cvfit$lambda.min)
# coefs = coefs[coefs != 0]
# control = coxph.control(iter.max = 0)
# multiCox=coxph(Surv(futime, fustat) ~ ., data = rt, control = control, init = coefs)
# 

multiCoxSum=summary(multiCox)

#输出模型参数
outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
write.table(outTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)

#输出病人风险值
riskScore=predict(multiCox, type="risk", newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime", "fustat", coxGene)

rt$risk <- riskScore
View(rt)


# Below use the Median as riskscore cutoff:
risk=as.vector(ifelse(riskScore>median(riskScore), "high", "low"))
outTab=cbind(rt[,outCol], riskScore=as.vector(riskScore), risk)
write.table(cbind(id=rownames(outTab),outTab), file="risk.txt", sep="\t", quote=F, row.names=F)

# Below use the best.optimal cutoff of the riskscore:
surv.cut_PRAD <- surv_cutpoint(
  rt,
  time = "futime",
  event = "fustat",
  variables = "risk"
)

summary(surv.cut_PRAD)
surv.cat_PRAD <- surv_categorize(surv.cut_PRAD)
View(surv.cat_PRAD)
diff <- survdiff(Surv(futime,fustat) ~ risk, data = surv.cat_PRAD)
pValue <- 1-pchisq(diff$chisq, df=1)
pValue <- signif(pValue, 4)
pValue <- format(pValue, scientific=FALSE)



################################################################################################################
################################# survival (immAS13.survival.R) ################################################
################################################################################################################


#引用包
library(survival)
library(survminer)
#setwd("C:\\biowolf\\immAS\\13.survival")     #设置工作目录

#定义生存曲线的函数
bioSurvival=function(inputFile=null, outFile=null){
  #读取输入文件
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  #比较高低风险组生存差异，得到显著性p值
  diff=survdiff(Surv(futime, fustat) ~ risk, data=rt)
  pValue=1-pchisq(diff$chisq, df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  #print(surv_median(fit))

  #绘制生存曲线
  surPlot=ggsurvplot(fit,
                     data=rt,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     surv.median.line = "hv",
                     legend.title="Risk",
                     legend.labs=c("High risk", "Low risk"),
                     xlab="Time(years)",
                     break.time.by = 1,
                     palette=c("red", "blue"),
                     risk.table=TRUE,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  pdf(file=outFile, onefile=FALSE, width=6.5, height=5.5)
  print(surPlot)
  dev.off()
}

#调用函数，绘制生存曲线
bioSurvival(inputFile="risk.txt", outFile="survival_PDCD2.pdf")




#Survival curve for MulCox with optimal riskscore cutoff----
#KM analysis and survival curve based on the high/low risk
fit_PRAD <- survfit(Surv(futime,fustat) ~ risk,
                    data = surv.cat_PRAD)

pdf(file = "survival_PDCD2.pdf", width = 8, height = 8)

surPlot_bestcut=ggsurvplot(fit_PRAD,
                           data=surv.cat_PRAD, #这个可要可不要
                           #conf.int=TRUE,
                           pval=TRUE, #显示log-rank test p-value
                           pval.size=6,
                           #legend.labs=c("High", "Low"), #最好别显示， 不然基因名会乱
                           #legend.title=gene,  #最好别显示， 不然基因名会乱
                           xlab="Time (years)",
                           ylab="Progression-free survival",
                           break.time.by = 1,
                           risk.table.title="Number at risk",
                           palette=c("#323232","#CCCCCC"),
                           risk.table=T,
                           risk.table.height=.25 )
surPlot_bestcut

dev.off()
#
#
# # riskScore table export for training set----
# ## re-define the risk to make sure to include corretly in the final 'risk' table
# risk <- surv.cat_PRAD$risk
# #risk <- as.vector(ifelse(riskScore>median(riskScore), "high","low")) #No more MEDIAN！！！
# View(rt)
# ##write out 'risk' table
# write.table(cbind(id = rownames(cbind(rt[ , outCol], riskScore, risk)), cbind(rt[ , outCol],riskScore,risk)),
#             file = "risk.txt",
#             sep = "\t",
#             quote = F,
#             row.names = F)

# # Covert the exported riskscore into linear format as inport the 'modfy' version----
# rt_GSig_PFS_PANCANCER_train <- read.table("/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Model construction/R_data_needed/risk_GSig_PFS_PANCANCER_345train_modfy.txt", header = T, sep = "\t", check.names = F, row.names = 1)
# #rt_GSig_PFS_PANCANCER_train <- read.table("/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/DFS.PFS/PFS_edgeR/split_kym/final_v7/risk_GSig_PFS_PANCANCER_308train.txt", header = T, sep = "\t", check.names = F, row.names = 1)
# View(rt_GSig_PFS_PANCANCER_train)




################################################################################################################
################################# riskPlot (immAS14.riskPlot.R) ################################################
################################################################################################################




library(pheatmap)       #引用包
#setwd("C:\\biowolf\\immAS\\14.riskPlot")      #设置工作目录

#定义风险曲线的函数
bioRiskPlot=function(inputFile=null, riskScoreFile=null, survStatFile=null, heatmapFile=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #读取输入文件
  rt$riskScore[rt$riskScore>quantile(rt$riskScore,0.99)]=quantile(rt$riskScore,0.99)
  rt$risk=factor(rt$risk, levels=c("low", "high"))
  rt=rt[order(rt$riskScore),]      #按照风险打分对样品排序
  
  #绘制风险曲线
  riskClass=rt[,"risk"]
  lowLength=length(riskClass[riskClass=="low"])
  highLength=length(riskClass[riskClass=="high"])
  lowMax=max(rt$riskScore[riskClass=="low"])
  line=rt[,"riskScore"]
  pdf(file=riskScoreFile, width=7, height=4)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing risk socre)", ylab="Risk score",
       col=c(rep("green",lowLength),rep("red",highLength)) )
  abline(h=lowMax,v=lowLength,lty=2)
  legend("topleft", c("High risk", "Low Risk"),bty="n",pch=19,col=c("red","green"),cex=1.2)
  dev.off()
  
  #绘制生存状态图
  color=as.vector(rt$fustat)
  color[color==1]="red"
  color[color==0]="green"
  pdf(file=survStatFile, width=7, height=4)
  plot(rt$futime, pch=19,
       xlab="Patients (increasing risk socre)", ylab="Survival time (years)",
       col=color)
  legend("topleft", c("Dead", "Alive"),bty="n",pch=19,col=c("red","green"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
  
  #绘制风险热图
  rt1=rt[c(3:(ncol(rt)-2))]
  rt1=t(rt1)
  annotation=data.frame(type=rt[,ncol(rt)])
  rownames(annotation)=rownames(rt)
  pdf(file=heatmapFile, width=7, height=4)
  pheatmap(rt1, 
           annotation=annotation, 
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           show_colnames = F,
           scale="row",
           color = colorRampPalette(c(rep("green",3), "white", rep("red",3)))(50),
           fontsize_col=3,
           fontsize=7,
           fontsize_row=8)
  dev.off()
}
#调用函数，绘制风险曲线
bioRiskPlot(inputFile="risk.txt",
            riskScoreFile="riskScore.pdf",
            survStatFile="survStat.pdf",
            heatmapFile="heatmap.pdf")







################################################################################################################
##################### independent signature assessment (immAS15.indep.R) #######################################
################################################################################################################




library(survival)       #引用包
#setwd("C:\\biowolf\\immAS\\15.indep")     #设置工作目录

############绘制森林图函数############
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
  #读取输入文件
  rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

  #输出图形
  pdf(file=forestFile, width=6.5, height=4.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))

  #绘制森林图左边的临床信息
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)

  #绘制右边的森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=3)
  abline(v=1, col="black", lty=2, lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=2)
  axis(1)
  dev.off()
}
############绘制森林图函数############

#定义独立预后分析函数
indep=function(riskFile=null,cliFile=null,uniOutFile=null,multiOutFile=null,uniForest=null,multiForest=null){
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #读取风险文件
  cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)      #读取临床文件

  #数据合并
  sameSample=intersect(row.names(cli),row.names(risk))
  risk=risk[sameSample,]
  cli=cli[sameSample,]
  rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])

  #单因素独立预后分析
  uniTab=data.frame()
  for(i in colnames(rt[,3:ncol(rt)])){
    cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
    coxSummary = summary(cox)
    uniTab=rbind(uniTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
  write.table(uniTab,file=uniOutFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=uniOutFile, forestFile=uniForest, forestCol="green")

  #多因素独立预后分析
  uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<1,]
  rt1=rt[,c("futime", "fustat", as.vector(uniTab[,"id"]))]
  multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
  multiCoxSum=summary(multiCox)
  multiTab=data.frame()
  multiTab=cbind(
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  multiTab=cbind(id=row.names(multiTab),multiTab)
  write.table(multiTab,file=multiOutFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=multiOutFile, forestFile=multiForest, forestCol="red")
}

#调用函数,进行独立预后分析; NOTE for TCGA-PRAD, there are three tumour samples duplicated with the same first three chunk id (xxxx-xx-xxxx), so need to identify and collapse somewhere else using the 'clinical.txt' in EXCEL, and then call the 'clinical.txt'; something like below:
# a <- read.table("clinical.txt", header=T, sep="\t", check.names=F)
# View(a)
# table(duplicated(a$id))
# b <- subset(a, duplicated(id))
# b
#
# # id Age Path_T Path_N Gleason
# # 54  TCGA-HC-8258  56     T2     N0       6
# # 74  TCGA-HC-7740  59     T2     N0       7
# # 236 TCGA-HC-8265  66     T3     N0       8

indep(riskFile="risk.txt",
      cliFile="clinical.txt",
      uniOutFile="uniCox_Clinical.txt",
      multiOutFile="multiCox_Clinical.txt",
      uniForest="uniForest.pdf",
      multiForest="multiForest.pdf")



################################################################################################################
############################################### ROC (immAS16.ROC.R) ############################################
################################################################################################################



#引用包
library(survival)
library(survminer)
library(timeROC)
riskFile="risk.txt"         #风险输入文件
cliFile="clinical.txt"      #临床数据文件
#setwd("C:\\biowolf\\immAS\\16.ROC")     #修改工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
View(risk)
risk=risk[,c("futime", "fustat", "riskScore")]

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
View(cli)
#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)
View(rt)

# NEW BIT ADDING FOR VARIABLE TRANSFORMATION (BUT DOES NOT MATTER FOR 2-LEVEL CATEGORICAL VARIABLE I THINK!!! ANSWER: YES!!!TRUE FOR BOTH UNI&MULTI)
# It needs to be number 1 or 0 or more or continuous (for sure) for ROC curve analysis
# rt[rt$age == '<60',]$age <- 0 
# rt[rt$age == '>=60',]$age <- 1 

rt[rt$Path_T == 'T2',]$Path_T <- 0 
rt[rt$Path_T == 'T3-4',]$Path_T <- 1 
rt$Path_T <-  as.numeric(rt$Path_T)


rt[rt$Path_N == 'N0',]$Path_N <- 0 
rt[rt$Path_N == 'N1',]$Path_N <- 1 
rt$Path_N<-  as.numeric(rt$Path_N)


rt[rt$Gleason == '<=7',]$Gleason <- 0 
rt[rt$Gleason == '>7',]$Gleason <- 1 
rt$Gleason <- as.numeric(rt$Gleason)


str(rt)

#定义颜色
bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)


######绘制1 3 5年的ROC曲线######
ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
               marker=risk$riskScore,cause=1,
               weighting='aalen',
               times=c(1,3,5),ROC=TRUE)
pdf(file="ROC_LASSOCOX_135yrs.pdf", width=5.5, height=5.5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=2)
plot(ROC_rt,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=bioCol[1:3], lwd=2, bty = 'n')
dev.off()


######绘制临床的ROC曲线######
predictTime=5    #定义预测年限
aucText=c()
pdf(file="cliROC_5yr.pdf", width=5.5, height=5.5)
#绘制风险得分的ROC曲线
i=3
ROC_rt=timeROC(T=risk$futime,
               delta=risk$fustat,
               marker=risk$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2)
aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
abline(0,1)
#对临床数据进行循环，绘制临床数据的ROC曲线
for(i in 4:ncol(rt)){
  ROC_rt=timeROC(T=rt$futime,
                 delta=rt$fustat,
                 marker=rt[,i], cause=1,
                 weighting='aalen',
                 times=c(predictTime),ROC=TRUE)
  plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2, add=TRUE)
  aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
}
#绘制图例，得到ROC曲线下的面积
legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])
dev.off()





################################################################################################################
############################################### Clinical Cor (immAS17.cliCor.R) ################################
################################################################################################################



#引用包
library(limma)
library(ggpubr)
riskFile="uniSigExp.txt"#"risk.txt"         #风险文件
cliFile="clinical_categorised.txt"      #临床数据文件
#setwd("C:\\biowolf\\immAS\\17.cliCor")     #设置工作目录

#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
View(risk)
risk$`CDC42BPA|10042|ES`[risk$`CDC42BPA|10042|ES`>quantile(risk$`CDC42BPA|10042|ES`,0.99)]=quantile(risk$`CDC42BPA|10042|ES`,0.99)

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk=risk[samSample,"CDC42BPA|10042|ES",drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk, cli)
View(rt)
#临床相关性分析，输出图形结果
for(clinical in colnames(rt[,2:ncol(rt)])){
  data=rt[c("CDC42BPA|10042|ES", clinical)]
  colnames(data)=c("CDC42BPA|10042|ES", "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  #设置比较组
  group=levels(factor(data$clinical))
  data$clinical=factor(data$clinical, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  #绘制箱线图
  boxplot=ggboxplot(data, x="clinical", y="CDC42BPA|10042|ES", color="clinical",
                    xlab=clinical,
                    ylab="PSI",
                    legend.title=clinical,
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons)
  #输出图片
  pdf(file=paste0(clinical, "_CDC42BPA|10042|ES.pdf"), width=5.5, height=5)
  print(boxplot)
  dev.off()
}






################################################################################################################
############################################### CNOMO and Cali (immAS18.Nomo.R) ################################
################################################################################################################






library(rms)                 #引用包
riskFile="risk.txt"          #风险输入文件
cliFile="clinical_categorised.txt"       #临床数据文件
#setwd("C:\\biowolf\\immAS\\18.Nomo")     #设置工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli=cli[apply(cli,1,function(x)any(is.na(match('unknow',x)))),,drop=F]

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)
paste(colnames(rt)[3:ncol(rt)],collapse="+")

#数据打包
dd <- datadist(rt)
options(datadist="dd")
#生成函数
#f <- cph(Surv(futime, fustat) ~ riskScore+Age+Gender+Grade+Stage+T+M+N, x=T, y=T, surv=T, data=rt, time.inc=1)
f <- cph(Surv(futime, fustat) ~ riskScore+Path_T+Path_N+Gleason, x=T, y=T, surv=T, data=rt, time.inc=1)
surv <- Survival(f)
#建立nomogram
nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(3, x), function(x) surv(5, x)), 
                lp=F, funlabel=c("1-year survival", "3-year survival", "5-year survival"), 
                maxscale=100, 
                fun.at=c(0.99, 0.9, 0.8, 0.7, 0.5, 0.3,0.1,0.01))  

#列线图可视化
pdf(file="Nomogram135_sigVariable.pdf", width=9.5, height=7.5)
plot(nom)
dev.off()

#校准曲线
time=5    #预测年限
f <- cph(Surv(futime, fustat) ~ riskScore+Age+Path_T+Path_N+Gleason, x=T, y=T, surv=T, data=rt, time.inc=time)
cal <- calibrate(f, cmethod="KM", method="boot", u=time, m=90, B=1000)
pdf(file="calibration_5y_sigVariable.pdf", width=9, height=8.5)
plot(cal,
     xlim=c(0,1),
     ylim=c(0,1),
     xlab=paste0("Nomogram-Predicted Probability of ", time, "-Year PFS"),
     ylab=paste0("Actual ", time, "-Year PFS(proportion)"), lwd=1.5,
     col="red", sub=T)
dev.off()



################################################################################################################
############################################### ROC of NOMO ####################################################
################################################################################################################



multiCox_nomo=coxph(Surv(futime, fustat) ~ riskScore+Path_T+Path_N+Gleason, data = rt)
multiCoxSum_nomo = summary(multiCox_nomo)
riskScore_nomo=predict(multiCox_nomo, type="risk", newdata=rt)
coxGene=rownames(multiCoxSum_nomo$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime", "fustat", coxGene)

rt$risk <- riskScore_nomo
View(rt)





# 
# #引用包
# library(survival)
# library(survminer)
# library(timeROC)
# riskFile="risk.txt"         #风险输入文件
# cliFile="clinical.txt"      #临床数据文件
# #setwd("C:\\biowolf\\immAS\\16.ROC")     #修改工作目录
# 
# #读取风险输入文件
# risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
# View(risk)
# risk=risk[,c("futime", "fustat", "riskScore")]
# 
# #读取临床数据文件
# cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
# View(cli)
# #合并数据
# samSample=intersect(row.names(risk), row.names(cli))
# risk1=risk[samSample,,drop=F]
# cli=cli[samSample,,drop=F]
# rt=cbind(risk1, cli)
# View(rt)

# NEW BIT ADDING FOR VARIABLE TRANSFORMATION (BUT DOES NOT MATTER FOR 2-LEVEL CATEGORICAL VARIABLE I THINK!!! ANSWER: YES!!!TRUE FOR BOTH UNI&MULTI)
# It needs to be number 1 or 0 or more or continuous (for sure) for ROC curve analysis
# rt[rt$age == '<60',]$age <- 0 
# rt[rt$age == '>=60',]$age <- 1 

rt[rt$Path_T == 'T2',]$Path_T <- 0 
rt[rt$Path_T == 'T3-4',]$Path_T <- 1 
rt$Path_T <-  as.numeric(rt$Path_T)


rt[rt$Path_N == 'N0',]$Path_N <- 0 
rt[rt$Path_N == 'N1',]$Path_N <- 1 
rt$Path_N<-  as.numeric(rt$Path_N)


rt[rt$Gleason == '<=7',]$Gleason <- 0 
rt[rt$Gleason == '>7',]$Gleason <- 1 
rt$Gleason <- as.numeric(rt$Gleason)


str(rt)

#定义颜色
bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)


######绘制1 3 5年的ROC曲线######
ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
               marker=rt$risk,cause=1,
               weighting='aalen',
               times=c(1,3,5),ROC=TRUE)
pdf(file="ROC_NOMOriskscore_135yrs.pdf", width=5.5, height=5.5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=2)
plot(ROC_rt,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=bioCol[1:3], lwd=2, bty = 'n')
dev.off()




