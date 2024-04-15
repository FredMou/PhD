
#================== this script uses the cli info directly from the PCaDB

setwd("/Users/zhuofanmou/Documents/PhD Exeter/ClariomD/alternative_splicing/R_workspace_v2/129immAS/workspace2/final_v3")
getwd()



TCGA_cli <- readRDS("/Users/zhuofanmou/Documents/PhD Exeter/ClariomD/alternative_splicing/R_workspace_v2/129immAS/workspace2/final/TCGA-PRAD_eSet.RDS")
TCGA_cli
TCGA_cli <- pData(TCGA_cli)
View(TCGA_cli)

TCGA_cli <- TCGA_cli[ , c("patient_id",
                          "sample_type",
                          "age_at_diagnosis",
                          "race",
                          "pathological_t_stage",
                          "pathological_n_stage",
                          "preop_psa",
                          "gleason_primary_pattern",
                          "gleason_secondary_pattern",
                          "gleason_score",
                          "time_to_death",
                          "os_status",
                          "time_to_bcr",
                          "bcr_status")]

TCGA_primary <- TCGA_cli[ TCGA_cli$sample_type == "Primary", ]
rownames(TCGA_primary) <- TCGA_primary$patient_id
View(TCGA_primary)

TCGA_survival <- TCGA_primary[ , c( "time_to_death",
                                    "os_status",
                                    "time_to_bcr",
                                    "bcr_status")]
TCGA_survival$time_to_death <- TCGA_survival$time_to_death/12
TCGA_survival$time_to_bcr <- TCGA_survival$time_to_bcr/12
View(TCGA_survival)





################################################################################
######### K-M Univariate Survival Analysis (log-rank test) for AS event ########
################################################################################
#Best cut-off point for continuous variables----
##Using surv_cutpoint:
##Determine the optimal cutpoint for one or multiple continuous variables at once, using the maximally selected rank statistics from the 'maxstat' R package. 
##This is an outcome-oriented methods providing a value of a cutpoint that correspond to the most significant relation with outcome (here, survival).
##Ref:
##https://www.jianshu.com/p/b5e96158059b?utm_campaign=hugo&utm_medium=reader_share&utm_content=note
##http://r-addict.com/2016/11/21/Optimal-Cutpoint-maxstat.html

###################### Library
library(survival)
library(survminer)

##################### Data import
# Event_dat=read.table( "asTime.txt", header=T, sep="\t", check.names=F, row.names=1)
Event_dat=load( "PSI_TCGA.rda")
Event_dat
Event_dat <- PSIexp
View(Event_dat)
Event_dat <- Event_dat[ , -1]

class(Event_dat)
dim(Event_dat)

Event_dat <- as.data.frame(Event_dat)
Event_dat$`PDCD2|78501|RI` <- as.numeric(Event_dat$`PDCD2|78501|RI`)
Event_dat$`ZWINT|11811|RI` <- as.numeric(Event_dat$`ZWINT|11811|RI` )
Event_dat$`TMUB2|41787|AT` <- as.numeric(Event_dat$`TMUB2|41787|AT`)
Event_dat$`PSMD8|49638|AT` <- as.numeric(Event_dat$`PSMD8|49638|AT`)
Event_dat$`ECE2|67860|AT` <- as.numeric(Event_dat$`ECE2|67860|AT`)
Event_dat$`CDC42BPA|10042|ES` <- as.numeric(Event_dat$`CDC42BPA|10042|ES`)
Event_dat$`TGIF1|44505|AT` <- as.numeric(Event_dat$`TGIF1|44505|AT`)
Event_dat$`EARS2|35597|RI` <- as.numeric(Event_dat$`EARS2|35597|RI` )
Event_dat$`SMAP1|76647|ES` <- as.numeric(Event_dat$`SMAP1|76647|ES`)
Event_dat$`TTLL5|28520|AT` <- as.numeric(Event_dat$`TTLL5|28520|AT`)
Event_dat$`MTMR14|63120|ES` <- as.numeric(Event_dat$`MTMR14|63120|ES`)
Event_dat$`SLC25A19|43430|AP` <- as.numeric(Event_dat$`SLC25A19|43430|AP`)
Event_dat$`TTC31|54078|RI`<- as.numeric(Event_dat$`TTC31|54078|RI`)
str(Event_dat)

##################### Rename the column names
colnames(Event_dat)[1:ncol(Event_dat)] <- c("PDCD2.78501.RI" , 
                              "ZWINT.11811.RI"  , 
                              "TMUB2.41787.AT"  , 
                              "PSMD8.49638.AT" ,  
                              "ECE2.67860.AT"   , 
                              "CDC42BPA.10042.ES",
                              "TGIF1.44505.AT" ,  
                              "EARS2.35597.RI"  ,  
                              "SMAP1.76647.ES" ,  
                              "TTLL5.28520.AT"   ,
                              "MTMR14.63120.ES" , 
                              "SLC25A19.43430.AP" ,
                              "TTC31.54078.RI" )




########################## Combine the data (inclding the PFS) with the BCR data
##
same=intersect(row.names(Event_dat),row.names(TCGA_survival)) 
length(same)

##

##
d1=Event_dat[same,] 
dim(d1)
head(rownames(d1))
d2 <- TCGA_survival[same,] 
dim(d2)
head(rownames(d2))
##


##Check if row names matched
table(row.names(d1) == row.names(d2))
##

## Combine the data (including PFS) for each of the cohorts
d3 <- cbind(d1,d2) 
str(d3)
d3 <- d3[ , c(14:17, 1:13)]
View(d3)

# BCR
table(d3$bcr_status)
d3_bcr <- d3[ , c(3,4,5:17)]
colnames(d3_bcr)[1] <- c("futime")
colnames(d3_bcr)[2] <- c("fustat")
View(d3_bcr)

# OS
table(d3$os_status)
d3_os <- d3[ , c(1,2,5:17)]
colnames(d3_os)[1] <- c("futime")
colnames(d3_os)[2] <- c("fustat")
View(d3_os)

######################### Output the made data for BCR analysis
TCGA_bcr <- d3_bcr
TCGA_bcr <- na.omit(TCGA_bcr)
sum(TCGA_bcr$futime<=0)
TCGA_bcr[TCGA_bcr$futime==0,] <- NA
TCGA_bcr <- na.omit(TCGA_bcr) # remove the rows with 0 survival time, so that LASSO-COX can run properly
View(TCGA_bcr)
table(TCGA_bcr$fustat)
save(TCGA_bcr, file ="TCGA_bcr.rda")

TCGA_os <- d3_os
TCGA_os <- na.omit(TCGA_os)
sum(TCGA_os$futime<=0)
TCGA_os[TCGA_os$futime==0,] <- NA
TCGA_os <- na.omit(TCGA_os) # remove the rows with 0 survival time, so that LASSO-COX can run properly
View(TCGA_os)
table(TCGA_os$fustat)
save(TCGA_os, file ="TCGA_os.rda")


#################################################################################################################
###################################### UniCox regression  ######################################################
################################################################################################################
#=========================================== For BCR
# library
library(survival)
#inputFile="asTime.txt"
pFilter=0.05

# read in data
#rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt <- TCGA_bcr
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
write.table(outTab,file="uniCox_limour_bcr.txt",sep="\t",row.names=F,quote=F)
#output for significant events from Uni Cox regression analysis
sigTab=outTab[as.numeric(as.vector(outTab$pvalue))<pFilter,]
write.table(sigTab,file="uniCox.Sig_limour_bcr.txt",sep="\t",row.names=F,quote=F)

##############Because we are now considering each event as a binary variable, so no need to use the PSI (continuous) values
# Uni-Cox significant AS events with their PSI
sigGenes=c("futime", "fustat", as.vector(sigTab[,1]))
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp), uniSigExp)
write.table(uniSigExp, file="uniSigExp_limour_bcr.txt", sep="\t", row.names=F, quote=F)


#=========================================== For OS
# library
library(survival)
#inputFile="asTime.txt"
pFilter=0.05

# read in data
#rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt <- TCGA_os
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
write.table(outTab,file="uniCox_limour_os.txt",sep="\t",row.names=F,quote=F)
#output for significant events from Uni Cox regression analysis
sigTab=outTab[as.numeric(as.vector(outTab$pvalue))<pFilter,]
write.table(sigTab,file="uniCox.Sig_limour_os.txt",sep="\t",row.names=F,quote=F)

##############Because we are now considering each event as a binary variable, so no need to use the PSI (continuous) values
# Uni-Cox significant AS events with their PSI
sigGenes=c("futime", "fustat", as.vector(sigTab[,1]))
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp), uniSigExp)
write.table(uniSigExp, file="uniSigExp_limour_os.txt", sep="\t", row.names=F, quote=F)



###############################################################################################################
########################################## LASSO-COX ##########################################################
###############################################################################################################
rt=read.table("uniSigExp_PCaDB_bcr.txt",header=T,sep="\t",check.names=F,row.names=1)


#rt <- Dat_BCR

#rt <- train

#rt <- Event_dat

View(rt)
library(glmnet)
# set a seed so lasso will be the same and reviewer can retrieve
set.seed(0)
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
coef <- coef(cvfit, s = cvfit$lambda.min)#coef(fit, s = cvfit$lambda.min) # use cvfit will be the same, why???
coef
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene=c("futime","fustat",lassoGene)
lassoSigExp=rt[,lassoGene]
lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
View(lassoSigExp)
write.table(lassoSigExp,file="lasso.sigExp.txt",sep="\t",row.names=F,quote=F)

######################################### LASSO-COX direct model prediction
# Ref: 
# https://github.com/rli012/PCaSignatures/blob/master/Evaluation_Survival_Analysis_Intra_Dataset.R
# https://github.com/rli012/PCaSignatures/blob/master/Evaluation_Survival_Analysis_Inter_Dataset.R
riskScore <- predict(cvfit, s=cvfit$lambda.min, newx = x, type="link") # response: reltive risk; link: linear prediction

outCol=lassoGene

rt$risk <- riskScore
View(rt)
class(rt)

library(maxstat)

opt_cut <-  maxstat.test(Surv(futime,fustat) ~ risk, data = rt, smethod = "LogRank",pmethod = "Lau92")#, pmethod = "Lau92")
opt_cut
cutpoint<- opt_cut$maxstats[[1]][4]
cutpoint

# # Below use the best cut-off as riskscore cutoff:
# training.risk.score <- riskScore[ , 1]
# class(training.risk.score)
# training.risk.group <- training.risk.score > cutpoint
# table(training.risk.group)


risk=as.vector(ifelse(riskScore>cutpoint, "high", "low"))
table(risk)
outTab=cbind(rt[,outCol], riskScore=as.vector(riskScore), risk)
write.table(cbind(id=rownames(outTab),outTab), file="risk.txt", sep="\t", quote=F, row.names=F)

risk=read.table("risk.txt", header=T, sep="\t", check.names=F, row.names=1)
View(risk)

diff <- survdiff(Surv(futime,fustat) ~ risk, data = risk)
pValue <- 1-pchisq(diff$chisq, df=1)
pValue <- signif(pValue, 4)
pValue <- format(pValue, scientific=FALSE)



fit_PRAD <- survfit(Surv(futime,fustat) ~ risk,
                    data = risk)

pdf(file = "survival_optimalBased_LASSOCOX.pdf", width = 8, height = 8)

surPlot_bestcut=ggsurvplot(fit_PRAD,
                           data=risk, #这个可要可不要
                           #conf.int=TRUE,
                           pval=TRUE, #显示log-rank test p-value
                           pval.size=6,
                           #legend.labs=c("High", "Low"), #最好别显示， 不然基因名会乱
                           #legend.title=gene,  #最好别显示， 不然基因名会乱
                           xlab="Time (years)",
                           ylab="Progression-free survival",
                           break.time.by = 1,
                           risk.table.title="Number at risk",
                           palette=c("#AB1815", "#27347F"),
                           risk.table=T,
                           risk.table.height=.25 )
surPlot_bestcut

dev.off()




########################################### ROC

#引用包
library(survival)
library(survminer)
library(timeROC)
riskFile="risk.txt"         #风险输入文件
cliFile="clinical_categorised.txt"      #临床数据文件
#setwd("C:\\biowolf\\immAS\\16.ROC")     #修改工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
View(risk)
risk=risk[,c("futime", "fustat", "riskScore")]

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
View(cli)
str(cli)
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

rt$Age<-ifelse(rt$Age==">=60",1,0)

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


##### time dependent ROC ######
ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
               marker=risk$riskScore,cause=1,
               weighting='aalen',
               times=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),ROC=TRUE)
ROC_rt$AUC


###### 1 3 5 8 10 year ROC######
pdf(file="ROC_LASSOCOX_135810yrs.pdf", width=5.5, height=5.5)
plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=2)
plot(ROC_rt,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=8,col=bioCol[4],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=10,col=bioCol[5],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[5])),
         paste0('AUC at 8 years: ',sprintf("%.03f",ROC_rt$AUC[8])),
         paste0('AUC at 10 years: ',sprintf("%.03f",ROC_rt$AUC[10]))),
       col=bioCol[1:5], lwd=2, bty = 'n')
dev.off()

###### 3 5 8 year ROC######
pdf(file="ROC_LASSOCOX_358yrs.pdf", width=5.5, height=5.5)
plot(ROC_rt,time=3,col=bioCol[1],title=FALSE,lwd=2)
plot(ROC_rt,time=5,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=8,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)

legend('bottomright',
       c(paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[5])),
         paste0('AUC at 8 years: ',sprintf("%.03f",ROC_rt$AUC[8]))),
       col=bioCol[1:3], lwd=2, bty = 'n')
dev.off()


######绘制临床的ROC曲线######
predictTime=10    #定义预测年限
aucText=c()
pdf(file="cliROC_10yr.pdf", width=5.5, height=5.5)
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



#c-index， similar to AUC， often used in COX evaluating cancer patients prognostic model's accuracy----
#BiocManager::install("survcomp")
library(survcomp)
# For training set
cindex <- concordance.index(x = risk$riskScore,
                            surv.time = risk$futime,
                            surv.event = risk$fustat,
                            method = "noether")

cindex$c.index
cindex$se
cindex$lower
cindex$upper
cindex$p.value

