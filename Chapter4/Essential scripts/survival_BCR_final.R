setwd("/Users/zhuofanmou/Downloads/survival_analysis")
getwd()


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
# # Event_dat=read.table( "asTime.txt", header=T, sep="\t", check.names=F, row.names=1)
# Event_dat=load( "PSI_TCGA.rda")
# Event_dat
# Event_dat <- PSIexp
# View(Event_dat)
# Event_dat <- Event_dat[ , -1]
# 
# class(Event_dat)
# dim(Event_dat)
# 
# Event_dat <- as.data.frame(Event_dat)
# Event_dat$`PDCD2|78501|RI` <- as.numeric(Event_dat$`PDCD2|78501|RI`)
# Event_dat$`ZWINT|11811|RI` <- as.numeric(Event_dat$`ZWINT|11811|RI` )
# Event_dat$`TMUB2|41787|AT` <- as.numeric(Event_dat$`TMUB2|41787|AT`)
# Event_dat$`PSMD8|49638|AT` <- as.numeric(Event_dat$`PSMD8|49638|AT`)
# Event_dat$`ECE2|67860|AT` <- as.numeric(Event_dat$`ECE2|67860|AT`)
# Event_dat$`CDC42BPA|10042|ES` <- as.numeric(Event_dat$`CDC42BPA|10042|ES`)
# Event_dat$`TGIF1|44505|AT` <- as.numeric(Event_dat$`TGIF1|44505|AT`)
# Event_dat$`EARS2|35597|RI` <- as.numeric(Event_dat$`EARS2|35597|RI` )
# Event_dat$`SMAP1|76647|ES` <- as.numeric(Event_dat$`SMAP1|76647|ES`)
# Event_dat$`TTLL5|28520|AT` <- as.numeric(Event_dat$`TTLL5|28520|AT`)
# Event_dat$`MTMR14|63120|ES` <- as.numeric(Event_dat$`MTMR14|63120|ES`)
# Event_dat$`SLC25A19|43430|AP` <- as.numeric(Event_dat$`SLC25A19|43430|AP`)
# Event_dat$`TTC31|54078|RI`<- as.numeric(Event_dat$`TTC31|54078|RI`)
# str(Event_dat)
# 
# 
# 
# clinical_BCR=read.csv("TCGA-PRAD_clinical_specific_modyfied.csv", row.names=1)
# rownames(clinical_BCR) <- clinical_BCR$bcr_patient_barcode
# clinical_BCR <- clinical_BCR[ , c(13,14,15,16)]
# clinical_BCR$BCR_time <- clinical_BCR$BCR_time/365
# clinical_BCR$os_time <- clinical_BCR$os_time/365
# table(clinical_BCR$BCR_status)
# table(clinical_BCR$os_status)
# View(clinical_BCR)
# str(clinical_BCR)
# 
# # # below for dcf limour
# # clinical_BCR <- clinical_BCR[ , -c(15,16)]
# # clinical_BCR$dcf_time <- clinical_BCR$dcf_time/365
# # clinical_BCR$os_time <- clinical_BCR$os_time/365
# # clinical_BCR <- clinical_BCR[ , c(11,12,13,14)]
# # table(clinical_BCR$dcf_status)
# # table(clinical_BCR$os_status)
# # View(clinical_BCR)
# # str(clinical_BCR)
# 
# 
# 
# ########################## Combine the data (inclding the PFS) with the BCR data
# ##
# same=intersect(row.names(Event_dat),row.names(clinical_BCR)) 
# length(same)
# 
# ##
# 
# ##
# d1=Event_dat[same,] 
# dim(d1)
# head(rownames(d1))
# d2 <- clinical_BCR[same,] 
# dim(d2)
# head(rownames(d2))
# ##
# 
# 
# ##Check if row names matched
# table(row.names(d1) == row.names(d2))
# ##
# 
# ## Combine the data (including PFS) for each of the cohorts
# d3 <- cbind(d1,d2) 
# str(d3)
# d3 <- d3[ , c(14:17,1:13)]
# View(d3)
# 
# 
# 
# 
# # #d3 <- d3[ , -c(3,4)]
# # colnames(d3)[1] <- c("futime")
# # colnames(d3)[2] <- c("fustat")
# # View(d3)
# # table(d3$fustat)
# # dim(d3)
# # class(d3)
# 
# ##################### Rename the column names
# colnames(d3)[5:ncol(d3)] <- c("PDCD2.78501.RI" , 
#                               "ZWINT.11811.RI"  , 
#                               "TMUB2.41787.AT"  , 
#                               "PSMD8.49638.AT" ,  
#                               "ECE2.67860.AT"   , 
#                               "CDC42BPA.10042.ES",
#                               "TGIF1.44505.AT" ,  
#                               "EARS2.35597.RI"  ,  
#                               "SMAP1.76647.ES" ,  
#                               "TTLL5.28520.AT"   ,
#                               "MTMR14.63120.ES" , 
#                               "SLC25A19.43430.AP" ,
#                               "TTC31.54078.RI" )
# 
# 
# 
# # BCR
# table(d3$BCR_status)
# d3_bcr <- d3[ , c(4,3,5:17)]
# colnames(d3_bcr)[1] <- c("futime")
# colnames(d3_bcr)[2] <- c("fustat")
# View(d3_bcr)
# 
# # OS
# table(d3$os_status)
# d3_os <- d3[ , c(2,1,5:17)]
# colnames(d3_os)[1] <- c("futime")
# colnames(d3_os)[2] <- c("fustat")
# View(d3_os)
# 
# ######################### Output the made data for BCR analysis
# TCGA_bcr <- d3_bcr
# TCGA_bcr <- na.omit(TCGA_bcr)
# sum(TCGA_bcr$futime<=0)
# TCGA_bcr[TCGA_bcr$futime==0,] <- NA
# TCGA_bcr <- na.omit(TCGA_bcr) # remove the rows with 0 survival time, so that LASSO-COX can run properly
# View(TCGA_bcr)    # !!! This is the same as "Dat_BCR.rda"
# table(TCGA_bcr$fustat)
# save(TCGA_bcr, file ="TCGA_bcr.rda")
# 
# TCGA_os <- d3_os
# TCGA_os <- na.omit(TCGA_os)
# sum(TCGA_os$futime<=0)
# TCGA_os[TCGA_os$futime==0,] <- NA
# TCGA_os <- na.omit(TCGA_os) # remove the rows with 0 survival time, so that LASSO-COX can run properly
# View(TCGA_os)
# table(TCGA_os$fustat)
# save(TCGA_os, file ="TCGA_os.rda")





######################### Output the made data for BCR analysis (OLD)
# Dat_BCR <- d3
# save(Dat_BCR, file ="Dat_BCR.rda")

########################################## This scrippt STARTS from HERE!!! #####################

##################################################################
############# Dataset split for training and testing #############
##################################################################



#==== read in BCRFS and AS event data
rt_BCRFS=read.table('asTime_BCRFS.txt', header=T, sep="\t", check.names=F, row.names=1)
View(rt_BCRFS)

#==== convert the column names to all '.', otherwise the output KM plot will not work properly!!!
colnames <- colnames(rt_BCRFS)
colnames
colnames_converted <- gsub("\\|", "\\.", colnames)
colnames_converted
colnames(rt_BCRFS) <- colnames_converted

#==== change the ONLY special gene name 'PPT2-EGFL8' to 'PPT2_EGFL8', other wise the output KM plot will not work properly!!!
which( colnames(rt_BCRFS)=="PPT2-EGFL8.75757.RI" ) # 91
colnames(rt_BCRFS)[91] <- 'PPT2_EGFL8.75757.RI' 
colnames(rt_BCRFS)[91]

#==== remove sample(s) than has survival time less than 1 month
rt_BCRFS <- rt_BCRFS[ -c(rt_BCRFS$futime < 1/12), ] # ONLY 1 sample identified and ,thus, removed

# # #==== remove one sample that not in the HTA dataset
# # View(rt_BCRFS)
# which( colnames(rt_BCRFS)=="CTNNB1.64249.RI" ) # 99
# rt_BCRFS <- rt_BCRFS[ , -99 ]
# 
# dim(rt_BCRFS)  # 412 samples * 142 columns (140 AS events)
# 


# Split data internally----
# https://www.jianshu.com/p/33425451f1f2?utm_campaign=hugo&utm_medium=reader_share&utm_content=note
library(caret)
set.seed(0109)
sam_BCRFS<- createDataPartition(rt_BCRFS$fustat, p = .7,list = FALSE)
train <- rt_BCRFS[sam_BCRFS,]
dim(train)
View(train)
#write.table(train, file = "trainingset.txt", sep = "\t", quote = F, row.names = F, col.names = F) # Export training set
save(train, file = "train.rda") # Export training set


test <- rt_BCRFS[-sam_BCRFS,]
dim(test)
#write.table(test, file = "testingset.txt", sep = "\t", quote = F, row.names = F, col.names = F) # Export teating set
save(test, file = "test.rda") # Export teating set





###############################################################################################
##################### KM and logrank test for each AS event (in batch) ########################
###############################################################################################

#==== KM and log-rank test analysis

# read in data
# For use as whole dataset
#train <- rt_BCRFS


# For use as splitted dataset
train <- train


View(train)
###################### best cut-ff approach for survival analysis
surv.cut_PRAD <- surv_cutpoint(
  train,
  time = "futime",
  event = "fustat",
  variables = c(colnames(train)[3:ncol(train)])
)


summary(surv.cut_PRAD)#cutpoint statistic 
#plot(surv.cut_PRAD, "MYC", palette = "npg") #MYC as an example

surv.cat_PRAD <- surv_categorize(surv.cut_PRAD) # Categorise each gene expression based on the optimal cut-off; if less=>low; if more=>high
View(surv.cat_PRAD)
class(surv.cat_PRAD)


###################### KM survival analysis and plot for single AS event (as batch)

picDir <- 'KM plot 141 AS events train dataset 0.7 vs 0.3_seed0109'
dir.create(picDir)
outTab <- data.frame()

#K-M plot for all 14 genes at once
for (gene in colnames(train[ ,3:ncol(train) ])){
  #a = train[ , gene] < median(train[ , gene])
  #a = CombindDat_model_fkpm[ , gene]<df_optcut[ 1, gene] # for optimal cut-off; 注意：这里就是前面所说的自己划分，没有用‘surv_categorize’ 公式， 同时也‘至少’多一个‘high’
  a=surv.cat_PRAD[ , gene] #这个是用‘surv_categorize’ 自带公式分high/low， 和上一行虽然少一个high但是得到的p值更加显著！且和之前得出来的五个基因一致！
  diff = survdiff(Surv(futime,fustat) ~ a, data = train)
  pValue = 1-pchisq(diff$chisq, df=1)
  outTab = rbind(outTab, cbind(gene=gene, pvalue=pValue))
  pValue=signif(pValue,4)
  pValue=format(pValue,scientific=FALSE)
  
  fit <- survfit(Surv(futime,fustat) ~ a, data =train)
  summary(fit)
  
  # tiff(file = paste(picDir, "/", gene, ".survival.tiff", sep = ""),
  #      width = 14,
  #      height = 14,
  #      units = "cm",
  #      compression = "lzw",
  #      bg = "white",
  #      res = 600)
  
  pdf(file = paste(picDir, "/", gene, ".survival.pdf", sep = ""), 
      width = 8, 
      height = 8)
  
  
  plot(fit,
       lwd = 2,
       col = c("#AB1815", "#27347F"),#c("#323232","#CCCCCC"),#col = c("red","blue"),
       xlab = "Time (year)",
       mark.time = T,
       ylab = "BCR-free survival"
       # main = paste("Progression-free survival curve (p=", pValue ,")", sep = "")
  )
  legend("topright",
         c(paste(gene," high expression",sep = ""),
           paste(gene," low expression",sep = "")),
         lwd = 2,
         col = c("#AB1815", "#27347F") #c("#323232","#CCCCCC")
  )
  legend("bottomright",
         paste("Logrank p=",pValue, sep = ""))
  dev.off()
}

##Export p-value table
write.table(outTab, file = "BCR_survival_km_141ASevents(single_optcut)_0.7_vs_0.3_seed0109.txt", sep = "\t", row.names = F, quote = F, col.names = T)


# ##If want, can chech the p-value manually; for example: HNRNPL
# diff = survdiff(Surv(futime,fustat) ~SRSF4.1424.ES, data = surv.cat_PRAD)
# pValue = 1-pchisq(diff$chisq, df=1)
# pValue=signif(pValue,4)
# pValue=format(pValue,scientific=TRUE)

SigKM_events <- outTab[outTab$pvalue<0.05 , ]
View(SigKM_events)
dim(SigKM_events)
rownames(SigKM_events) <- SigKM_events$gene

#################################################################################################################
###################################### UniCox regression  ######################################################
################################################################################################################

pFilter=0.05

# read in data
# For use as whole dataset
#rt <- rt_BCRFS


# For use as splitted dataset
rt <- train


str(rt)
class(rt)
View(rt)
#==== Univariate Cox regression in batch
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
write.table(outTab,file="uniCox_limour_BCRonly.txt",sep="\t",row.names=F,quote=F)

#output for significant events from Uni Cox regression analysis
sigTab=outTab[as.numeric(as.vector(outTab$pvalue))<pFilter,]
write.table(sigTab,file="uniCox.Sig_limour_BCRonly.txt",sep="\t",row.names=F,quote=F)

# Uni-Cox significant AS events with their PSI
sigGenes=c("futime", "fustat", as.vector(sigTab[,1]))
uniSigExp=rt[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp), uniSigExp)
write.table(uniSigExp, file="uniSigExp_limour_BCRonly.txt", sep="\t", row.names=F, quote=F)

###############################################################################################################
########################################## Taking the intersection ##########################################################
###############################################################################################################

View(SigKM_events)
UniSig_events <- uniSigExp
View(UniSig_events)
class(UniSig_events)
UniSig_events <- UniSig_events[ , -c(1)]

#Combine2
same=intersect(row.names(SigKM_events),colnames(UniSig_events))
length(same)


d1=UniSig_events[,c(1,2)]
View(d1)
dim(d1)
head(rownames(d1))
d2 <- UniSig_events[,same]
View(d2)
dim(d2)
head(rownames(d2))

##Check if row names matched
table(row.names(d1) == row.names(d2))

d3 <- cbind(d1,d2)
combined <- d3
# combined$patient_id <- rownames(combined)
View(combined)
write.table(combined,"SigEvent_intersect_UniCox_and_KM.txt",quote=F,row.names = F,col.names = T,sep="\t")







###############################################################################################################
########################################## LASSO-COX ##########################################################
###############################################################################################################

## For not using KM and UniCox intersected set
#rt=read.table("uniSigExp_limour_BCRonly.txt",header=T,sep="\t",check.names=F,row.names=1)

## For using KM and UniCox interested set
rt <- combined

View(rt)
library(glmnet)
# set a seed so lasso will be the same and reviewer can retrieve
set.seed(0109)
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

# ######################################### LASSO-COX direct model prediction
# # Ref: 
# # https://github.com/rli012/PCaSignatures/blob/master/Evaluation_Survival_Analysis_Intra_Dataset.R
# # https://github.com/rli012/PCaSignatures/blob/master/Evaluation_Survival_Analysis_Inter_Dataset.R
# riskScore <- predict(cvfit, s=cvfit$lambda.min, newx = x, type="link") # response: reltive risk; link: linear prediction
# 
# outCol=lassoGene
# 
# rt$risk <- riskScore
# View(rt)
# class(rt)
# 
# library(maxstat)
# 
# opt_cut <-  maxstat.test(Surv(futime,fustat) ~ risk, data = rt, smethod = "LogRank",pmethod = "Lau92")#, pmethod = "Lau92")
# opt_cut
# cutpoint<- opt_cut$maxstats[[1]][4]
# cutpoint
# 
# # # Below use the best cut-off as riskscore cutoff:
# # training.risk.score <- riskScore[ , 1]
# # class(training.risk.score)
# # training.risk.group <- training.risk.score > cutpoint
# # table(training.risk.group)
# 
# 
# risk=as.vector(ifelse(riskScore>cutpoint, "high", "low"))
# table(risk)
# outTab=cbind(rt[,outCol], riskScore=as.vector(riskScore), risk)
# write.table(cbind(id=rownames(outTab),outTab), file="risk.txt", sep="\t", quote=F, row.names=F)
# 
# risk=read.table("risk.txt", header=T, sep="\t", check.names=F, row.names=1)
# View(risk)
# 
# diff <- survdiff(Surv(futime,fustat) ~ risk, data = risk)
# pValue <- 1-pchisq(diff$chisq, df=1)
# pValue <- signif(pValue, 4)
# pValue <- format(pValue, scientific=FALSE)
# 
# 
# 
# fit_PRAD <- survfit(Surv(futime,fustat) ~ risk,
#                     data = risk)
# 
# pdf(file = "survival_optimalBased_LASSOCOX_limour_BCRonly.pdf", width = 8, height = 8)
# 
# surPlot_bestcut=ggsurvplot(fit_PRAD,
#                            data=risk, #这个可要可不要
#                            #conf.int=TRUE,
#                            pval=TRUE, #显示log-rank test p-value
#                            pval.size=6,
#                            #legend.labs=c("High", "Low"), #最好别显示， 不然基因名会乱
#                            #legend.title=gene,  #最好别显示， 不然基因名会乱
#                            xlab="Time (years)",
#                            ylab="Progression-free survival",
#                            break.time.by = 1,
#                            risk.table.title="Number at risk",
#                            palette=c("#AB1815", "#27347F"),
#                            risk.table=T,
#                            risk.table.height=.25 )
# surPlot_bestcut
# 
# dev.off()
# 
# 
# 
# 
# ########################################### ROC
# 
# #引用包
# library(survival)
# library(survminer)
# library(timeROC)
# riskFile="risk.txt"         #风险输入文件
# cliFile="clinical_categorised.txt"      #临床数据文件
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
# str(cli)
# #合并数据
# samSample=intersect(row.names(risk), row.names(cli))
# risk1=risk[samSample,,drop=F]
# cli=cli[samSample,,drop=F]
# rt=cbind(risk1, cli)
# View(rt)
# 
# # NEW BIT ADDING FOR VARIABLE TRANSFORMATION (BUT DOES NOT MATTER FOR 2-LEVEL CATEGORICAL VARIABLE I THINK!!! ANSWER: YES!!!TRUE FOR BOTH UNI&MULTI)
# # It needs to be number 1 or 0 or more or continuous (for sure) for ROC curve analysis
# # rt[rt$age == '<60',]$age <- 0
# # rt[rt$age == '>=60',]$age <- 1
# 
# rt$Age<-ifelse(rt$Age==">=60",1,0)
# 
# rt[rt$Path_T == 'T2',]$Path_T <- 0
# rt[rt$Path_T == 'T3-4',]$Path_T <- 1
# rt$Path_T <-  as.numeric(rt$Path_T)
# 
# 
# rt[rt$Path_N == 'N0',]$Path_N <- 0
# rt[rt$Path_N == 'N1',]$Path_N <- 1
# rt$Path_N<-  as.numeric(rt$Path_N)
# 
# 
# rt[rt$Gleason == '<=7',]$Gleason <- 0
# rt[rt$Gleason == '>7',]$Gleason <- 1
# rt$Gleason <- as.numeric(rt$Gleason)
# 
# 
# str(rt)
# 
# #定义颜色
# bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)
# 
# 
# ##### time dependent ROC ######
# ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
#                marker=risk$riskScore,cause=1,
#                weighting='aalen',
#                times=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),ROC=TRUE)
# ROC_rt$AUC
# 
# 
# ###### 1 3 5 8 10 year ROC######
# pdf(file="ROC_LASSOCOX_135810yrs_limour_BCRonly.pdf", width=5.5, height=5.5)
# plot(ROC_rt,time=1,col=bioCol[1],title=FALSE,lwd=2)
# plot(ROC_rt,time=3,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
# plot(ROC_rt,time=5,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
# plot(ROC_rt,time=8,col=bioCol[4],add=TRUE,title=FALSE,lwd=2)
# plot(ROC_rt,time=10,col=bioCol[5],add=TRUE,title=FALSE,lwd=2)
# legend('bottomright',
#        c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
#          paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3])),
#          paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[5])),
#          paste0('AUC at 8 years: ',sprintf("%.03f",ROC_rt$AUC[8])),
#          paste0('AUC at 10 years: ',sprintf("%.03f",ROC_rt$AUC[10]))),
#        col=bioCol[1:5], lwd=2, bty = 'n')
# dev.off()
# 
# ###### 3 5 8 year ROC######
# pdf(file="ROC_LASSOCOX_358yrs_limour_BCRonly.pdf", width=5.5, height=5.5)
# plot(ROC_rt,time=3,col=bioCol[1],title=FALSE,lwd=2)
# plot(ROC_rt,time=5,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
# plot(ROC_rt,time=8,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
# 
# legend('bottomright',
#        c(paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3])),
#          paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[5])),
#          paste0('AUC at 8 years: ',sprintf("%.03f",ROC_rt$AUC[8]))),
#        col=bioCol[1:3], lwd=2, bty = 'n')
# dev.off()
# 
# 
# ######绘制临床的ROC曲线######
# predictTime=8    #定义预测年限
# aucText=c()
# pdf(file="cliROC_8yr_limoue_BCRonly.pdf", width=5.5, height=5.5)
# #绘制风险得分的ROC曲线
# i=3
# ROC_rt=timeROC(T=risk$futime,
#                delta=risk$fustat,
#                marker=risk$riskScore, cause=1,
#                weighting='aalen',
#                times=c(predictTime),ROC=TRUE)
# plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2)
# aucText=c(paste0("Risk", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
# abline(0,1)
# #对临床数据进行循环，绘制临床数据的ROC曲线
# for(i in 4:ncol(rt)){
#   ROC_rt=timeROC(T=rt$futime,
#                  delta=rt$fustat,
#                  marker=rt[,i], cause=1,
#                  weighting='aalen',
#                  times=c(predictTime),ROC=TRUE)
#   plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2, add=TRUE)
#   aucText=c(aucText, paste0(colnames(rt)[i],", AUC=",sprintf("%.3f",ROC_rt$AUC[2])))
# }
# #绘制图例，得到ROC曲线下的面积
# legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(rt)-1)])
# dev.off()
# 
# 
# 
# #c-index， similar to AUC， often used in COX evaluating cancer patients prognostic model's accuracy----
# #BiocManager::install("survcomp")
# library(survcomp)
# # For training set
# cindex <- concordance.index(x = risk$riskScore,
#                             surv.time = risk$futime,
#                             surv.event = risk$fustat,
#                             method = "noether")
# 
# cindex$c.index
# cindex$se
# cindex$lower
# cindex$upper
# cindex$p.value




################################################################################################################
########################################### Multivariate-COX ###################################################
################################################################################################################


#==== importing LASSO screened events
#rt <- combined
rt <- read.table("lasso.sigExp.txt",header=T,sep="\t",check.names=F,row.names=1) # use the events selected by LASSO
#rt=read.table("uniSigExp_limour_BCRonly.txt",header=T,sep="\t",check.names=F,row.names=1) # use the events selected by univariate COX
View(rt)

multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)#coxph(Surv(futime, fustat) ~ `PDCD2|78501|RI`, data = rt)#coxph(Surv(futime, fustat) ~ ., data = rt)


# ## stepwise regression: https://www.jianshu.com/p/33425451f1f2?utm_campaign=hugo&utm_medium=reader_share&utm_content=note
#
# # colnames(rt) <- c("futime","fustat","TTC31.54078.RI",
# #               "PSMD8.49638.AT",
# #               "PDCD2.78501.RI",
# #               "CDC42BPA.10042.ES",
# #               "TGIF1.44505.AT")
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
multiCoxSum

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

#==== riskscore
# REF: https://stats.stackexchange.com/questions/209030/predict-from-a-cox-model-with-beta-coefficients

# #== use relative risk
# riskScore=predict(multiCox, type="risk", newdata=rt)
# View(as.data.frame(riskScore))

#== use linear prediction
# linear formulation for COX
## Ref: https://github.com/rli012/PCaSignatures/blob/master/Evaluation_Survival_Analysis_Intra_Dataset.R

## For training set
riskScore<- predict(multiCox, type="lp", newdata=rt, se.fit = FALSE)

## For testing set
rt <- test
riskScore<- predict(multiCox, type="lp", newdata=rt, se.fit = FALSE)

## For whole set
rt <- rt_BCRFS
riskScore<- predict(multiCox, type="lp", newdata=rt, se.fit = FALSE)


# IMPORTANT, in many paper, they use the traditional beta*X to present riskscore, the the riskScore just caluculated need to be converted a little bit to match the same expression
which( colnames(rt)=="CYP4F12.48110.RI" | colnames(rt)=="NFATC4.26991.RI" | colnames(rt)=="PIGO.86233.RI" | colnames(rt)=="CYP3A5.80711.RI" | colnames(rt)=="ALS2CL.64461.RI"  | colnames(rt)=="FXYD3.49039.RI")
X <- rt[ , c(3,5,7,8,10,11)] # For training
X <- rt[ , c(12,26,72,101,116,133)] # For testing and complete
View(X)
riskScore <- riskScore + as.vector(multiCox$coefficients %*% as.matrix(colMeans(X)))
class(riskScore)
View(as.data.frame(riskScore))

#== use mannually Bets*X, linear prediction: NOTE! this is different to the LP above as above uses Beta*(X-mean(X)):https://stats.stackexchange.com/questions/209030/predict-from-a-cox-model-with-beta-coefficients

# riskScore_mannual <- 1.293*(rt$CYP4F12.48110.RI-mean(rt$CYP4F12.48110.RI)) + 2.09*(rt$ZNF154.52329.RI-mean(rt$ZNF154.52329.RI)) + 0.9199*(rt$CYP3A5.80711.RI-mean(rt$CYP3A5.80711.RI)) + 1.206*(rt$ALS2CL.64461.RI-mean(rt$ALS2CL.64461.RI)) + 11.69*(rt$FXYD3.49039.RI-mean(rt$FXYD3.49039.RI))  # This is the SAME as 'predict(multiCox, type="lp", newdata=rt, se.fit = FALSE)'
# riskScore_mannual <- 1.293*(rt$CYP4F12.48110.RI) + 2.09*(rt$ZNF154.52329.RI) + 0.9199*(rt$CYP3A5.80711.RI) + 1.206*(rt$ALS2CL.64461.RI) + 11.69*(rt$FXYD3.49039.RI) # This is the SAME as 'as.matrix(X) %*% multiCox$coefficients'
# riskScore_mannual <- as.matrix(X) %*% multiCox$coefficients
# View(as.data.frame(riskScore_mannual))



coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime", "fustat", coxGene)

rt$risk <- riskScore
View(rt)


# # Below use the Median as riskscore cutoff:
# risk=as.vector(ifelse(riskScore>median(riskScore), "high", "low"))
# outTab=cbind(rt[,outCol], riskScore=as.vector(riskScore), risk)
# write.table(cbind(id=rownames(outTab),outTab), file="risk.txt", sep="\t", quote=F, row.names=F)

# Below use the best.optimal cutoff of the riskscore:
surv.cut_PRAD <- surv_cutpoint(
  rt,
  time = "futime",
  event = "fustat",
  variables = "risk"
)

summary(surv.cut_PRAD)
#maxstat.test(Surv(futime,fustat) ~ risk, data = rt, smethod = "LogRank", pmethod = "Lau92") # alternative check
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

# #定义生存曲线的函数; use this codes if used median as cut-off above
# bioSurvival=function(inputFile=null, outFile=null){
#   #读取输入文件
#   rt=read.table(inputFile, header=T, sep="\t", check.names=F)
#   #比较高低风险组生存差异，得到显著性p值
#   diff=survdiff(Surv(futime, fustat) ~ risk, data=rt)
#   pValue=1-pchisq(diff$chisq, df=1)
#   if(pValue<0.001){
#     pValue="p<0.001"
#   }else{
#     pValue=paste0("p=",sprintf("%.03f",pValue))
#   }
#   fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
#   #print(surv_median(fit))
#
#   #绘制生存曲线
#   surPlot=ggsurvplot(fit,
#                      data=rt,
#                      conf.int=T,
#                      pval=pValue,
#                      pval.size=6,
#                      surv.median.line = "hv",
#                      legend.title="Risk",
#                      legend.labs=c("High risk", "Low risk"),
#                      xlab="Time(years)",
#                      break.time.by = 1,
#                      palette=c("red", "blue"),
#                      risk.table=TRUE,
#                      risk.table.title="",
#                      risk.table.col = "strata",
#                      risk.table.height=.25)
#   pdf(file=outFile, onefile=FALSE, width=6.5, height=5.5)
#   print(surPlot)
#   dev.off()
# }
#
# #调用函数，绘制生存曲线
# bioSurvival(inputFile="risk.txt", outFile="survival_medianBased_CAT.pdf")




#Survival curve for MulCox with optimal riskscore cutoff----
#KM analysis and survival curve based on the high/low risk
fit_PRAD <- survfit(Surv(futime,fustat) ~ risk,
                    data = surv.cat_PRAD)

pdf(file = "survival_optimalBased_UCOXLASSOMCOX_limour_BCRonly_train.pdf", width = 8, height = 8)

surPlot_bestcut=ggsurvplot(fit_PRAD,
                           data=surv.cat_PRAD, #这个可要可不要
                           #conf.int=TRUE,
                           pval=TRUE, #显示log-rank test p-value
                           pval.size=6,
                           #legend.labs=c("High", "Low"), #最好别显示， 不然基因名会乱
                           #legend.title=gene,  #最好别显示， 不然基因名会乱
                           xlab="Time (years)",
                           ylab="BCR-free survival",
                           break.time.by = 1,
                           risk.table.title="Number at risk",
                           palette=c("#AB1815", "#27347F"),
                           risk.table=T,
                           risk.table.height=.25 )
surPlot_bestcut

dev.off()




# riskScore table export for training set----
## re-define the risk to make sure to include corretly in the final 'risk' table
risk <- surv.cat_PRAD$risk
#risk <- as.vector(ifelse(riskScore>median(riskScore), "high","low")) #No more MEDIAN！！！
#View(risk)
##write out 'risk' table
write.table(cbind(id = rownames(cbind(rt[ , outCol], riskScore, risk)), cbind(rt[ , outCol],riskScore,risk)),
            file = "risk.txt",
            sep = "\t",
            quote = F,
            row.names = F)






############################## ROC
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


######绘制 N 年的ROC曲线######
ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
               marker=risk$riskScore,cause=1,
               weighting='aalen',
               times=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),ROC=TRUE)
ROC_rt$AUC
# ####### For 3 and 5-year ROC only
# pdf(file="ROC_5810yrs.pdf", width=5.5, height=5.5)
# plot(ROC_rt,time=5,col="#27347F",title=FALSE,lwd=2)
# plot(ROC_rt,time=8,col="#AB1815",add=TRUE,title=FALSE,lwd=2)
# plot(ROC_rt,time=10,col="#008000",add=TRUE,title=FALSE,lwd=2)
# legend('bottomright',
#        c(paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[5])),
#          paste0('AUC at 8 years: ',sprintf("%.03f",ROC_rt$AUC[8])),
#          paste0('AUC at 10 years: ',sprintf("%.03f",ROC_rt$AUC[10]))),
# 
#        col=c("#27347F","#AB1815","#008000"), lwd=2, bty = 'n')
# dev.off()


###### 1 3 5 8 10 year ROC######
pdf(file="ROC_LASSOCOX_MULCOX_135810yrs_limour_BCRonly.pdf", width=5.5, height=5.5)
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
pdf(file="ROC_LASSOCOX_MULCOX_358yrs_limour_BCRonly_complete.pdf", width=5.5, height=5.5)
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
predictTime=8   #定义预测年限
aucText=c()
pdf(file="ROC_8yr_SignatureCli_limour_BCRonly.pdf", width=5.5, height=5.5)
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

# > cindex$c.index      
# [1] 0.6617601
# > cindex$se
# [1] 0.03595452
# > cindex$lower
# [1] 0.5912905
# > cindex$upper
# [1] 0.7322296
# > cindex$p.value
# [1] 6.826774e-06


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
      cliFile="clinical_categorised.txt",
      uniOutFile="uniCox_Clinical.txt",
      multiOutFile="multiCox_Clinical.txt",
      uniForest="uniForest.pdf",
      multiForest="multiForest.pdf")



################################################################################################################
############################################### Clinical Cor (immAS17.cliCor.R) ################################
################################################################################################################



#引用包
library(limma)
library(ggpubr)
riskFile="risk.txt"         #风险文件
cliFile="clinical_categorised.txt"      #临床数据文件
#setwd("C:\\biowolf\\immAS\\17.cliCor")     #设置工作目录

# #读取风险文件
# risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
# View(risk)
# risk$`CDC42BPA|10042|ES`[risk$`CDC42BPA|10042|ES`>quantile(risk$`CDC42BPA|10042|ES`,0.99)]=quantile(risk$`CDC42BPA|10042|ES`,0.99)
# 
# #读取临床数据文件
# cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
# 
# #合并数据
# samSample=intersect(row.names(risk), row.names(cli))
# risk=risk[samSample,"CDC42BPA|10042|ES",drop=F]
# cli=cli[samSample,,drop=F]
# rt=cbind(risk, cli)
# View(rt)
# #临床相关性分析，输出图形结果
# for(clinical in colnames(rt[,2:ncol(rt)])){
#   data=rt[c("CDC42BPA|10042|ES", clinical)]
#   colnames(data)=c("CDC42BPA|10042|ES", "clinical")
#   data=data[(data[,"clinical"]!="unknow"),]
#   #设置比较组
#   group=levels(factor(data$clinical))
#   data$clinical=factor(data$clinical, levels=group)
#   comp=combn(group,2)
#   my_comparisons=list()
#   for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
#   #绘制箱线图
#   boxplot=ggboxplot(data, x="clinical", y="CDC42BPA|10042|ES", color="clinical",
#                     xlab=clinical,
#                     ylab="PSI",
#                     legend.title=clinical,
#                     add = "jitter")+ 
#     stat_compare_means(comparisons = my_comparisons)
#   #输出图片
#   pdf(file=paste0(clinical, "_CDC42BPA|10042|ES.pdf"), width=5.5, height=5)
#   print(boxplot)
#   dev.off()
# }




#读取风险文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk$riskScore[risk$riskScore>quantile(risk$riskScore,0.99)]=quantile(risk$riskScore,0.99)

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk=risk[samSample,"riskScore",drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk, cli)

#临床相关性分析，输出图形结果
for(clinical in colnames(rt[,2:ncol(rt)])){
  data=rt[c("riskScore", clinical)]
  colnames(data)=c("riskScore", "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  #设置比较组
  group=levels(factor(data$clinical))
  data$clinical=factor(data$clinical, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  #绘制箱线图
  boxplot=ggboxplot(data, x="clinical", y="riskScore", color="clinical",
                    xlab=clinical,
                    ylab="Risk score",
                    legend.title=clinical,
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons)
  #输出图片
  pdf(file=paste0(clinical, ".pdf"), width=5.5, height=5)
  print(boxplot)
  dev.off()
}





################################################################################################################
################################################ NOMO and Cali (immAS18.Nomo.R) ################################
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
#f <- cph(Surv(futime, fustat) ~ riskScore+Path_T+Path_N+Gleason, x=T, y=T, surv=T, data=rt, time.inc=1)
f <- cph(Surv(futime, fustat) ~ riskScore+Path_T+Gleason, x=T, y=T, surv=T, data=rt, time.inc=1)

surv <- Survival(f)
#建立nomogram
nom <- nomogram(f, fun=list(function(x) surv(3, x), function(x) surv(5, x), function(x) surv(8, x)), 
                lp=F, funlabel=c("3-year survival", "5-year survival", "8-year survival"), 
                maxscale=100, 
                fun.at=c(0.99, 0.9, 0.8, 0.7, 0.5, 0.3,0.1,0.01))  

#列线图可视化
pdf(file="Nomogram358_sigVariable_LASSOMCOX_limour_BCRonly.pdf", width=9.5, height=7.5)
plot(nom)
dev.off()

#校准曲线
time=8    #预测年限
f <- cph(Surv(futime, fustat) ~ riskScore+Path_T+Gleason, x=T, y=T, surv=T, data=rt, time.inc=time)
cal <- calibrate(f, cmethod="KM", method="boot", u=time, m=115, B=1000)
pdf(file="calibration_8y_sigVariable_LASSOMCOX_limour_BCRonly.pdf", width=9, height=8.5)
plot(cal,
     xlim=c(0,1),
     ylim=c(0,1),
     xlab=paste0("Nomogram-Predicted Probability of ", time, "-Year BCR"),
     ylab=paste0("Actual ", time, "-Year BCR(proportion)"), lwd=1.5,
     col="red", sub=T)
dev.off()

################################################################################################################
############################################### ROC of NOMO ####################################################
################################################################################################################

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

multiCox_nomo=coxph(Surv(futime, fustat) ~ riskScore+Path_T+Gleason, data = rt)

multiCoxSum_nomo = summary(multiCox_nomo)
riskScore_nomo=predict(multiCox_nomo, type="lp", newdata=rt)

# IMPORTANT, in many paper, they use the traditional beta*X to present riskscore, the the riskScore just caluculated need to be converted a little bit to match the same expression
which( colnames(rt)=="riskScore" | colnames(rt)=="Path_T" | colnames(rt)=="Gleason")
X <- rt[ , c(3,5,7)] 

View(X)
riskScore_nomo <- riskScore_nomo + as.vector(multiCox_nomo$coefficients %*% as.matrix(colMeans(X)))
class(riskScore_nomo)
View(as.data.frame(riskScore_nomo))



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

# rt[rt$Path_T == 'T2',]$Path_T <- 0
# rt[rt$Path_T == 'T3-4',]$Path_T <- 1
# rt$Path_T <-  as.numeric(rt$Path_T)
#
#
# rt[rt$Path_N == 'N0',]$Path_N <- 0
# rt[rt$Path_N == 'N1',]$Path_N <- 1
# rt$Path_N<-  as.numeric(rt$Path_N)
#
#
# rt[rt$Gleason == '<=7',]$Gleason <- 0
# rt[rt$Gleason == '>7',]$Gleason <- 1
# rt$Gleason <- as.numeric(rt$Gleason)
#
#
# str(rt)

#定义颜色
bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)


######绘制1 3 5年的ROC曲线######
ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
               marker=rt$risk,cause=1,
               weighting='aalen',
               times=c(3,5,8),ROC=TRUE)
ROC_rt$AUC
pdf(file="ROC_NOMOriskscore_358yrs_limour_BCRonly.pdf", width=5.5, height=5.5)
plot(ROC_rt,time=3,col=bioCol[1],title=FALSE,lwd=2)
plot(ROC_rt,time=5,col=bioCol[2],add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=8,col=bioCol[3],add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 8 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=bioCol[1:3], lwd=2, bty = 'n')
dev.off()


######绘制临床的ROC曲线######
predictTime=5   #定义预测年限
aucText=c()
pdf(file="ROC_3yr_NomoCli_limour_BCRonly.pdf", width=5.5, height=5.5)
#绘制风险得分的ROC曲线
i=3
ROC_rt=timeROC(T=rt$futime,
               delta=rt$fustat,
               marker=rt$riskScore, cause=1,
               weighting='aalen',
               times=c(predictTime),ROC=TRUE)
plot(ROC_rt, time=predictTime, col=bioCol[i-2], title=FALSE, lwd=2)
aucText=c(paste0("Event signature", ", AUC=", sprintf("%.3f",ROC_rt$AUC[2])))
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
# cindex <- concordance.index(x = risk$riskScore,
#                             surv.time = risk$futime,
#                             surv.event = risk$fustat,
#                             method = "noether")

cindex <- concordance.index(x = rt$risk,
                            surv.time = rt$futime,
                            surv.event = rt$fustat,
                            method = "noether")

cindex$c.index
cindex$se
cindex$lower
cindex$upper
cindex$p.value



# > cindex$c.index
# [1] 0.7237771
# > cindex$se
# [1] 0.03143967
# > cindex$lower
# [1] 0.6621565
# > cindex$upper
# [1] 0.7853977
# > cindex$p.value
# [1] 1.097679e-12

























