setwd('/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Model construction/R_data_needed')
getwd()


##################################################################
########################## Library needed ########################
##################################################################
library(limma) 
library(dplyr)
library(survival)
library(survminer)
library(glmnet)
library(ROCR)
library(caret)
library(Hmisc) # for c-index
library(timeROC)
library(ggplot2)
library(survcomp)
library(pheatmap)
library(pROC)
library(wesanderson)
library(stringr)
library(rms)
library(survivalROC)
##################################################################
########################## Loading the dataset ###################
##################################################################
lnames_pfs_pan_edgeR <- load("/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Model construction/R_data_needed/TCGA_PAN_PFS_GSig_Input_edgeR.rda") #load('/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/DFS.PFS/PFS_edgeR/split_kym/final_v4/TCGA_PAN_PFS_GSig_Input_edgeR.rda')
lnames_pfs_pan_edgeR
CombindDat_model_pfs_edgeR <-rt_final_PFS_edgeR
#CombindDat_model_pfs_deseq2$months <- CombindDat_model_pfs_deseq2$days_to_last_follow_up*365/30
View(CombindDat_model_pfs_edgeR)
CombindDat_model_pfs_edgeR$vital_status <- as.numeric(CombindDat_model_pfs_edgeR$vital_status)
str(CombindDat_model_pfs_edgeR)
table(CombindDat_model_pfs_edgeR$vital_status)

CombindDat_model_pfs_edgeR <- CombindDat_model_pfs_edgeR[CombindDat_model_pfs_edgeR$days_to_last_follow_up>=1/12,] #0.25
dim(CombindDat_model_pfs_edgeR)
View(CombindDat_model_pfs_edgeR)
##################################################################
############# Dataset split for training and testing #############
##################################################################

# Split data internally----
# https://www.jianshu.com/p/33425451f1f2?utm_campaign=hugo&utm_medium=reader_share&utm_content=note
# library(caret)
set.seed(123)
sam_pfs<- createDataPartition(CombindDat_model_pfs_edgeR$vital_status, p = .7,list = FALSE)
train <- CombindDat_model_pfs_edgeR[sam_pfs,]
dim(train)
View(train)
#write.table(train, file = "trainingset.txt", sep = "\t", quote = F, row.names = F, col.names = F) # Export training set
save(train, file = "train.rda") # Export training set


test <- CombindDat_model_pfs_edgeR[-sam_pfs,]
dim(test)
#write.table(test, file = "testingset.txt", sep = "\t", quote = F, row.names = F, col.names = F) # Export teating set
save(test, file = "test.rda") # Export teating set




##################################################################################
############################## Univariate Cox ####################################
### (Only one gene significant; maybe too stringent; Not necessary from Jack) ####
##################################################################################
#Uni Cox----
outTab_UniCOX <- data.frame()

table(train$vital_status)

for(i in colnames(train[, 3:ncol(train)])){
  cox <- coxph(Surv(days_to_last_follow_up,vital_status) ~ train[ , i], data = train)
  coxSummary = summary(cox)
  outTab_UniCOX <- rbind(outTab_UniCOX, cbind(gene=i, HR=coxSummary$coefficients[,"exp(coef)"], 
                                              z=coxSummary$coefficients[,"z"], 
                                              pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
}

write.table(outTab_UniCOX, file = "UnivariateCox_10genes_PFS_PAN_edgeR_345trainset.txt", sep = "\t", row.names = F, quote = F,col.names = T) 




##################################################################
######### K-M Survival Curve (log-rank test) for Gene Exp ########
##################################################################
#Best cut-off point for continuous variables----
##Using surv_cutpoint:
##Determine the optimal cutpoint for one or multiple continuous variables at once, using the maximally selected rank statistics from the 'maxstat' R package. 
##This is an outcome-oriented methods providing a value of a cutpoint that correspond to the most significant relation with outcome (here, survival).
##Ref:
##https://www.jianshu.com/p/b5e96158059b?utm_campaign=hugo&utm_medium=reader_share&utm_content=note
##http://r-addict.com/2016/11/21/Optimal-Cutpoint-maxstat.html


# best cut-ff approach for survival analysis

surv.cut_PRAD <- surv_cutpoint(
  train,
  time = "days_to_last_follow_up",
  event = "vital_status",
  variables = c('HNRNPL', 'ELAVL1', 'PCBP1', 'PCBP2', 'PABPN1','PTPRF','MRPL24', 'DANCR', 'MYC', 'TRPM4')
)


summary(surv.cut_PRAD)#cutpoint statistic MYC 5.637017  2.802201
#plot(surv.cut_PRAD, "MYC", palette = "npg") #MYC as an example

surv.cat_PRAD <- surv_categorize(surv.cut_PRAD) # Categorise each gene expression based on the optimal cut-off; if less=>low; if more=>high
View(surv.cat_PRAD)

#KM survival analysis and plot for single gene (as batch)----

picDir <- 'KM plot 10genes edgeR (opt-cut) final PFS from TCGA PANCANCER 345trainset'
dir.create(picDir)


outTab <- data.frame()

#K-M plot for all 14 genes at once
for (gene in colnames(train[ ,3:ncol(train) ])){
  #a = train[ , gene] < median(train[ , gene])
  #a = CombindDat_model_fkpm[ , gene]<df_optcut[ 1, gene] # for optimal cut-off; 注意：这里就是前面所说的自己划分，没有用‘surv_categorize’ 公式， 同时也‘至少’多一个‘high’
  a=surv.cat_PRAD[ , gene] #这个是用‘surv_categorize’ 自带公式分high/low， 和上一行虽然少一个high但是得到的p值更加显著！且和之前得出来的五个基因一致！
  diff = survdiff(Surv(days_to_last_follow_up,vital_status) ~ a, data = train)
  pValue = 1-pchisq(diff$chisq, df=1)
  outTab = rbind(outTab, cbind(gene=gene, pvalue=pValue))
  pValue=signif(pValue,4)
  pValue=format(pValue,scientific=FALSE)
  
  fit <- survfit(Surv(days_to_last_follow_up,vital_status) ~ a, data =train)
  summary(fit)
  
  tiff(file = paste(picDir, "/", gene, ".survival.tiff", sep = ""),
       width = 14,
       height = 14,
       units = "cm",
       compression = "lzw",
       bg = "white",
       res = 600)
  
  plot(fit,
       lwd = 2,
       col = c("#323232","#CCCCCC"),#col = c("red","blue"),
       xlab = "Time (year)",
       mark.time = T,
       ylab = "Progression-free survival"
       # main = paste("Progression-free survival curve (p=", pValue ,")", sep = "")
  )
  legend("topright",
         c(paste(gene," high expression",sep = ""),
           paste(gene," low expression",sep = "")),
         lwd = 2,
         col = c("#323232","#CCCCCC"))
  legend("bottomright",
         paste("Logrank p=",pValue, sep = ""))
  dev.off()
}

##Export p-value table
write.table(outTab, file = "survival_km_10genes(single_optcut)_edgeR_final_PFS_PANCANCER_345trainset.txt", sep = "\t", row.names = F, quote = F, col.names = T)


##If want, can chech the p-value manually; for example: HNRNPL
# diff = survdiff(Surv(days_to_last_follow_up,vital_status) ~HNRNPL, data = surv.cat_PRAD)
# pValue = 1-pchisq(diff$chisq, df=1)
# pValue=signif(pValue,4)
# pValue=format(pValue,scientific=TRUE)



#Survival curve for MulCox Option2 (with risk table)----
#KM analysis and survival curve based on the high/low risk
# fit_PRAD <- survfit(Surv(days_to_last_follow_up,vital_status) ~ PABPN1,
#                     data = surv.cat_PRAD)
# 
# pdf(file = "survival_MulCox_KMplot_GSig_edgeR_pfs_pan_308train_withRiskTab.pdf", width = 8, height = 8)
# 
# surPlot_bestcut=ggsurvplot(fit_PRAD,
#                            data=surv.cat_PRAD, #这个可要可不要
#                            #conf.int=TRUE,
#                            pval=TRUE, #显示log-rank test p-value
#                            pval.size=6,
#                            #legend.labs=c("High", "Low"), #最好别显示， 不然基因名会乱
#                            #legend.title=gene,  #最好别显示， 不然基因名会乱
#                            xlab="Time (years)",
#                            ylab="Progression-free survival",
#                            break.time.by = 1,
#                            risk.table.title="Number at risk",
#                            palette=c("#323232","#CCCCCC"),
#                            risk.table=T,
#                            risk.table.height=.25 )
# surPlot_bestcut
# 
# dev.off()


##################################################################
############## LASSO regression for gene selection ###############
############ (Use this as first selection step by Jack) ##########
##################################################################
#Ref:
#1.https://www.jianshu.com/p/33425451f1f2?utm_campaign=hugo&utm_medium=reader_share&utm_content=note
#2.KYM course codes

#Selecting the optimal/minimal Lambda (only! without split to train&test&ROC)----
#Select the genes of interest (PROBABLY from KM log-rank test above), here the five genes
x_lasso_edgeR_train <- as.matrix(train[ ,c('HNRNPL', 'ELAVL1', 'PCBP1', 'PCBP2', 'PABPN1','PTPRF','MRPL24', 'DANCR', 'MYC', 'TRPM4')])#as.matrix(train[ ,3:ncol(train) ])
dim(x_lasso_edgeR_train)
View(x_lasso_edgeR_train)
class(x_lasso_edgeR_train)
#Select the survial status and time for lasso survival (SXZXW):
y_lasso_edgeR_train<- data.matrix(Surv(train$days_to_last_follow_up, train$vital_status))
class(y_lasso_edgeR_train)
length(y_lasso_edgeR_train)

# library(glmnet)
set.seed(123)
#LASSO (KYM)
fitCV <- cv.glmnet(x=x_lasso_edgeR_train, y=y_lasso_edgeR_train,alpha=1, family = "cox",maxit=1000) # train x and train y' , maxit=1000
############## Plot1: for lambda ##########
# tiff(file = "lasso_lambda1.tiff",
#      width = 20,
#      height = 12,
#      units = "cm",
#      compression = "lzw",
#      bg = "white",
#      res = 600)

#PDF
pdf(file = "lasso_lambda1_edgeR_PFS_pan.pdf", width = 8, height = 8)
plot(fitCV)
abline(v=log(c(fitCV$lambda.min, fitCV$lambda.1se)), lty = "dashed")
dev.off()
####
#Lambda
lambda.min=fitCV$lambda.min
lambda.min
#Get the coefficient
coef.min = coef(fitCV, s = "lambda.min")
coef.min
#Check the coefficient
fit <- glmnet(x=x_lasso_edgeR_train, y=y_lasso_edgeR_train, family = "cox", alpha = 1,maxit = 1000) # , lambda = fitCV$lambda.min; make sure alpha = 1; 'lambda = lasso.train$lambda.min' is NOT necessary

############ Plot2: for coefficient#####
# tiff(file = "lasso_lambda2.tiff",
#      width = 20,
#      height = 12,
#      units = "cm",
#      compression = "lzw",
#      bg = "white",
#      res = 600)

#PDF
pdf(file = "lasso_lambda2_edgeR_PFS_pan.pdf", width = 8, height = 8)

plot(fit, xvar="lambda",label = TRUE)
# abline(v = log(lambda.min,10), lty = 3,
#        lwd = 2,
#        col = "black")
dev.off()
####



coefficients = as.matrix(coef(fit,fitCV$lambda.min))
idx = which(coefficients != 0)
coefficients[idx,]

actCoef <- coefficients[idx]
lassoGene <- rownames(coefficients)[idx]
write.table(lassoGene, file = "lassoGene.txt", sep = "\t", quote = F, row.names = F, col.names = F) # Export list of genes screened from LASSO



##################################################################
####################### MulCox construction ######################
##################################################################
# Inputting the selected gene signature from LASSO----
cox_Mul <- coxph(Surv(days_to_last_follow_up,vital_status) ~ PTPRF   +   DANCR    +    MYC     + TRPM4   +   PCBP1   +  PABPN1 , data = train)
cox_Mul <- coxph(Surv(days_to_last_follow_up,vital_status) ~  PCBP1+PABPN1+PTPRF+DANCR+MYC+TRPM4, data = train)
# Didirection stepwise screening step----
cox_Mul <- step(cox_Mul, direction = "both") #TRPM4 removed
summary <- summary(cox_Mul)
coxGene <- rownames(summary$coefficients)
coxGene <- gsub("`","",coxGene)
outCol <- c("days_to_last_follow_up","vital_status",coxGene)


# Risk scores (hazard functions) for training and testing sets----
riskScore_Gsig_PFS_PANCANCER_train <- predict(cox_Mul, type = "risk", newdata = train)# we need to use this riskscore （continuous）later on for classifying high&low risk groups
riskScore_Gsig_PFS_PANCANCER_test <- predict(cox_Mul, type = "risk", newdata = test)# we need to use this riskscore （continuous）later on for classifying high&low risk groups


#check C-index----
# library(Hmisc)
1-with(train, rcorr.cens(riskScore_Gsig_PFS_PANCANCER_train, Surv(days_to_last_follow_up,vital_status)))[1] #0.681
1-with(test, rcorr.cens(riskScore_Gsig_PFS_PANCANCER_test, Surv(days_to_last_follow_up,vital_status)))[1] #0.660


# optiml cut-off to classify riskscore on the Train and Test set now!----
train$risk <- riskScore_Gsig_PFS_PANCANCER_train
View(train)
test$risk <- riskScore_Gsig_PFS_PANCANCER_test
View(test)


# Optimal cut-ff approach for survival analysis of calculated riskScore (TRAINING cohort)----
surv.cut_PRAD <- surv_cutpoint(
  train,
  time = "days_to_last_follow_up",
  event = "vital_status",
  variables = "risk"
)

summary(surv.cut_PRAD)
surv.cat_PRAD <- surv_categorize(surv.cut_PRAD)
View(surv.cat_PRAD)
diff <- survdiff(Surv(days_to_last_follow_up,vital_status) ~ risk, data = surv.cat_PRAD)
pValue <- 1-pchisq(diff$chisq, df=1)
pValue <- signif(pValue, 4)
pValue <- format(pValue, scientific=FALSE)


#Plot: Overall/progression-free survival for GSig----
fit <- survfit(Surv(days_to_last_follow_up,vital_status) ~ risk, data = surv.cat_PRAD)
summary(fit) #check for 5-year survival rate

# tiff(file = "survival_MulCox_KMplot_MYC&ELAVL1_fpkm.tiff",
#      width = 14,
#      height = 14,
#      units = "cm",
#      compression = "lzw",
#      bg = "white",
#      res = 600)
#
#PDF
# pdf(file = "survival_MulCox_KMplot_GSig_edgeR_pfs_pan_308train.pdf", width = 8, height = 8)
# 
plot(fit,
     lwd = 2,
     col = c("#323232","#CCCCCC"),#col = c("red","blue"),
     xlab = "Time (year)",
     mark.time = T,
     ylab = "Progression-free survival"
     # main = paste("Progression-free survival curve (p=", pValue ,")", sep = "")
)
legend("topright",
       c("high risk", "low risk"),
       lwd = 2,
       col = c("#323232","#CCCCCC"))
legend("bottomright",
       paste("Logrank p=",pValue, sep = ""))
# 
# dev.off()
####


#Survival curve for MulCox Option2 (with risk table)----
#KM analysis and survival curve based on the high/low risk
fit_PRAD <- survfit(Surv(days_to_last_follow_up,vital_status) ~ risk,
                    data = surv.cat_PRAD)

pdf(file = "survival_MulCox_KMplot_GSig_edgeR_pfs_pan_345train_withRiskTab.pdf", width = 8, height = 8)

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


# riskScore table export for training set----
## re-define the risk to make sure to include corretly in the final 'risk' table
risk <- surv.cat_PRAD$risk
#risk <- as.vector(ifelse(riskScore>median(riskScore), "high","low")) #No more MEDIAN！！！
View(train)
##write out 'risk' table
write.table(cbind(id = rownames(cbind(train[ , outCol], riskScore_Gsig_PFS_PANCANCER_train, risk)), cbind(train[ , outCol],riskScore_Gsig_PFS_PANCANCER_train,risk)),
            file = "risk_GSig_PFS_PANCANCER_345train.txt",
            sep = "\t",
            quote = F,
            row.names = F)

# MulCox results export for training set----
write.table(cbind(id=coxGene, summary$coefficients),
            file = "MultivariateCox_coef_GSig_pfs_pan_345TRAIN.txt", 
            sep = "\t",
            quote = F,
            row.names = F)
cox_Mul


# Covert the exported riskscore into linear format as inport the 'modfy' version----
rt_GSig_PFS_PANCANCER_train <- read.table("/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Model construction/R_data_needed/risk_GSig_PFS_PANCANCER_345train_modfy.txt", header = T, sep = "\t", check.names = F, row.names = 1)
#rt_GSig_PFS_PANCANCER_train <- read.table("/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/DFS.PFS/PFS_edgeR/split_kym/final_v7/risk_GSig_PFS_PANCANCER_308train.txt", header = T, sep = "\t", check.names = F, row.names = 1)
View(rt_GSig_PFS_PANCANCER_train)

# MulCox forest----
# tiff(file = "forest_MulCox_GSig_deseq2_pfs_pan_218TRAIN.tiff",
#      width = 20,
#      height = 12,
#      units = "cm",
#      compression = "lzw",
#      bg = "white",
#      res = 600)
#PDF
pdf(file = "forest_MulCox_GSig_edgeR_pfs_pan_345TRAIN.pdf", width = 12, height = 8)

ggforest(cox_Mul,
         data=train,
         main = "Hazard ratio of PRAD",
         cpositions = c(0.02, 0.22,0.4),
         fontsize = 0.7,
         refLabel = "reference",
         noDigits = 2)

dev.off()



# Optimal cut-ff approach for survival analysis of calculated riskScore (TESTING cohort)----

surv.cut_PRAD <- surv_cutpoint(
  test,
  time = "days_to_last_follow_up",
  event = "vital_status",
  variables = "risk"
)

summary(surv.cut_PRAD)
surv.cat_PRAD <- surv_categorize(surv.cut_PRAD)
View(surv.cat_PRAD)
diff <- survdiff(Surv(days_to_last_follow_up,vital_status) ~ risk, data = surv.cat_PRAD)
pValue <- 1-pchisq(diff$chisq, df=1)
pValue <- signif(pValue, 4)
pValue <- format(pValue, scientific=FALSE)



# Plot: Overall/progression-free survival for GSig----
fit <- survfit(Surv(days_to_last_follow_up,vital_status) ~ risk, data = surv.cat_PRAD)
summary(fit) #check for 5-year survival rate
# tiff(file = "survival_MulCox_KMplot_MYC&ELAVL1_fpkm.tiff",
#      width = 14,
#      height = 14,
#      units = "cm",
#      compression = "lzw",
#      bg = "white",
#      res = 600)
#
#PDF
# pdf(file = "survival_MulCox_KMplot_GSig_edgeR_pfs_pan_165test.pdf", width = 8, height = 8)
# 
# 
# 
plot(fit,
     lwd = 2,
     col = c("#323232","#CCCCCC"),#col = c("red","blue"),
     xlab = "Time (year)",
     mark.time = T,
     ylab = "Progression-free survival"
     # main = paste("Progression-free survival curve (p=", pValue ,")", sep = "")
)
legend("topright",
       c("high risk", "low risk"),
       lwd = 2,
       col = c("#323232","#CCCCCC"))
legend("bottomright",
       paste("Logrank p=",pValue, sep = ""))
# 
# dev.off()
####


# Survival curve for MulCox Option2 (with risk table)----
#KM analysis and survival curve based on the high/low risk
fit_PRAD <- survfit(Surv(days_to_last_follow_up,vital_status) ~ risk,
                    data = surv.cat_PRAD)


pdf(file = "survival_MulCox_KMplot_GSig_edgeR_pfs_pan_147test_withRiskTab.pdf", width = 8, height = 8)

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

#riskScore table export for testing set----
## re-define the risk to make sure to include corretly in the final 'risk' table
risk <- surv.cat_PRAD$risk
#risk <- as.vector(ifelse(riskScore>median(riskScore), "high","low")) #No more MEDIAN！！！
View(test)
##write out 'risk' table
write.table(cbind(id = rownames(cbind(test[ , outCol], riskScore_Gsig_PFS_PANCANCER_test, risk)), cbind(test[ , outCol],riskScore_Gsig_PFS_PANCANCER_test,risk)),
            file = "risk_GSig_PFS_PANCANCER_147test.txt",
            sep = "\t",
            quote = F,
            row.names = F)

# Covert the exported riskscore into linear format as inport the 'modfy' version----
#rt_GSig_PFS_PANCANCER_test <- read.table("/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/DFS.PFS/PFS_edgeR/split_kym/final_v7/risk_GSig_PFS_PANCANCER_165test.txt", header = T, sep = "\t", check.names = F, row.names = 1)

rt_GSig_PFS_PANCANCER_test <- read.table("/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Model construction/R_data_needed/risk_GSig_PFS_PANCANCER_147test_mody.txt", header = T, sep = "\t", check.names = F, row.names = 1)
View(rt_GSig_PFS_PANCANCER_test)


##################################################################
################## Time-dependent ROC analysis ###################
##################################################################
#From: https://www.jianshu.com/p/b5e96158059b?utm_campaign=hugo&utm_medium=reader_share&utm_content=note
## For training set----
View(rt_GSig_PFS_PANCANCER_train)
colnames(rt_GSig_PFS_PANCANCER_train)
result <-with(rt_GSig_PFS_PANCANCER_train, timeROC( T=days_to_last_follow_up,
                                                    delta=vital_status,
                                                    marker=riskScore_Gsig_PFS_PANCANCER_train,#marker=riskScore_Gsig_PFS_PANCANCER_train_modfy,
                                                    cause=1,
                                                    times=c(1,3,5),
                                                    iid = TRUE))
identical(c(result$TP[,1],result$TP[,2],result$TP[,3]),as.numeric(result$TP))
dat = data.frame(fpr = as.numeric(result$FP),
                 tpr = as.numeric(result$TP),
                 time = rep(as.factor(c(1,3,5)),each = nrow(result$TP)))



# #TIFF
# tiff(file = "MulCox_ROC_allyrs_2Genesig_fpkm.tiff",
#      width = 20,
#      height = 20,
#      units = "cm",
#      compression = "lzw",
#      bg = "white",
#      res = 600)

#PDF for training set
pdf(file = "ROC_135yrs_Genesig_345TCGAvalid.pdf", width = 8, height = 8)


par(oma = c(0.5,1,0,1), font.lab = 1.5, font.axis = 1.5)

ggplot() +
  geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) +
  scale_color_manual(name = NULL,values = c("#323232", "#989898", "#CCCCCC"),
                     labels = paste0("AUC of ",c(1,3,5),"-year survival: ",
                                     format(round(result$AUC,3),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()+
  ggtitle("ROC curve of gene signature for PFS in TCGA-PRAD training cohort")+
  theme(plot.title = element_text(hjust = 0.5,colour = "black",face = "bold.italic",size = 10))

dev.off()

## Alt ROC2
# library(survivalROC)
# 
# 
# pdf(file="ROC.pdf",width=6,height=6)
# par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
# roc=survivalROC(Stime=rt_GSig_PFS_PANCANCER_train$days_to_last_follow_up, status=rt_GSig_PFS_PANCANCER_train$vital_status, marker = rt_GSig_PFS_PANCANCER_train$riskScore_Gsig_PFS_PANCANCER_train,
#                 predict.time =8, method="KM")
# plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red',
#      xlab="False positive rate", ylab="True positive rate",
#      main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
#      lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
# abline(0,1)
# dev.off()




## For testing set----
View(rt_GSig_PFS_PANCANCER_test)
colnames(rt_GSig_PFS_PANCANCER_test)
result <-with(rt_GSig_PFS_PANCANCER_test, timeROC( T=days_to_last_follow_up,
                                                   delta=vital_status,
                                                   marker=riskScore_Gsig_PFS_PANCANCER_test,#marker=riskScore_Gsig_PFS_PANCANCER_test_modfy,
                                                   cause=1,
                                                   times=c(1,3,5),
                                                   iid = TRUE))
identical(c(result$TP[,1],result$TP[,2],result$TP[,3]),as.numeric(result$TP))
dat = data.frame(fpr = as.numeric(result$FP),
                 tpr = as.numeric(result$TP),
                 time = rep(as.factor(c(1,3,5)),each = nrow(result$TP)))

#PDF for testing set
pdf(file = "ROC_135yrs_Genesig_147TCGAvalid.pdf", width = 8, height = 8)


par(oma = c(0.5,1,0,1), font.lab = 1.5, font.axis = 1.5)

ggplot() +
  geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) +
  scale_color_manual(name = NULL,values = c("#323232", "#989898", "#CCCCCC"),
                     labels = paste0("AUC of ",c(1,3,5),"-year survival: ",
                                     format(round(result$AUC,3),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()+
  ggtitle("ROC curve of gene signature for PFS in TCGA-PRAD testing cohort")+
  theme(plot.title = element_text(hjust = 0.5,colour = "black",face = "bold.italic",size = 10))

dev.off()
##

# Alt ROC2
# library(survivalROC)
# 
# 
# pdf(file="ROC.pdf",width=6,height=6)
# par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
# roc=survivalROC(Stime=rt_GSig_PFS_PANCANCER_test$days_to_last_follow_up, status=rt_GSig_PFS_PANCANCER_test$vital_status, marker = rt_GSig_PFS_PANCANCER_test$riskScore_Gsig_PFS_PANCANCER_test,
#                 predict.time =3, method="KM") #predict 3 year
# plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red',
#      xlab="False positive rate", ylab="True positive rate",
#      main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
#      lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
# abline(0,1)
# dev.off()


#c-index， similar to AUC， often used in COX evaluating cancer patients prognostic model's accuracy----
#BiocManager::install("survcomp")
library(survcomp)
# For training set
cindex <- concordance.index(x = rt_GSig_PFS_PANCANCER_train$riskScore_Gsig_PFS_PANCANCER_train,
                            surv.time = rt_GSig_PFS_PANCANCER_train$days_to_last_follow_up,
                            surv.event = rt_GSig_PFS_PANCANCER_train$vital_status,
                            method = "noether")

cindex$c.index
cindex$se
cindex$lower
cindex$upper
cindex$p.value

#For testing set
cindex <- concordance.index(x = rt_GSig_PFS_PANCANCER_test$riskScore_Gsig_PFS_PANCANCER_test,
                            surv.time = rt_GSig_PFS_PANCANCER_test$days_to_last_follow_up,
                            surv.event = rt_GSig_PFS_PANCANCER_test$vital_status,
                            method = "noether")

cindex$c.index
cindex$se
cindex$lower
cindex$upper
cindex$p.value



##################################################################
################ Calibrate curve for GSig Only ###################
##################################################################
## Ref: https://blog.csdn.net/zhongkeyuanchongqing/article/details/118622251?utm_medium=distribute.pc_relevant_t0.none-task-blog-2%7Edefault%7ECTRLIST%7Edefault-1.no_search_link&depth_1-utm_source=distribute.pc_relevant_t0.none-task-blog-2%7Edefault%7ECTRLIST%7Edefault-1.no_search_link

## For training cohort----
View(rt_GSig_PFS_PANCANCER_train)


ddlist <- datadist(rt_GSig_PFS_PANCANCER_train)
options(datadist = 'ddlist')
units(rt_GSig_PFS_PANCANCER_train$days_to_last_follow_up) <- "year"


pdf(file = "Cali_GSig_1yr_345train.pdf", width = 13, height = 8)
# 1-year calibration
f1_train <- cph(Surv(days_to_last_follow_up, vital_status) ~ PCBP1+PABPN1+PTPRF+DANCR+MYC, x=T, y=T, surv=T, data=rt_GSig_PFS_PANCANCER_train, time.inc = 1)
cal1_train <- calibrate(f1_train, cmethod = 'KM', method = 'boot', u=1,m=115, B=1000)
plot(cal1_train, xlab="Predicted PFS 1 year", ylab="Observed PFS 1 year")
dev.off()

pdf(file = "Cali_GSig_3yr_345train.pdf", width = 13, height = 8)
# 3-year calibration
f3_train <- cph(Surv(days_to_last_follow_up, vital_status) ~ PCBP1+PABPN1+PTPRF+DANCR+MYC, x=T, y=T, surv=T, data=rt_GSig_PFS_PANCANCER_train, time.inc = 3) # riskScore_Gsig_PFS_PANCANCER_train
cal3_train <- calibrate(f3_train, cmethod = 'KM', method = 'boot', u=3,m=115, B=1000)
plot(cal3_train,xlab="Predicted PFS 3 years", ylab="Observed PFS 3 years")
dev.off()

pdf(file = "Cali_GSig_5yr_345train.pdf", width = 13, height = 8)
# 5-year calibration
f5_train <- cph(Surv(days_to_last_follow_up, vital_status) ~ PCBP1+PABPN1+PTPRF+DANCR+MYC, x=T, y=T, surv=T, data=rt_GSig_PFS_PANCANCER_train, time.inc = 5)
cal5_train <- calibrate(f5_train, cmethod = 'KM', method = 'boot', u=5,m=115, B=1000)
plot(cal5_train,xlab="Predicted PFS 5 years", ylab="Observed PFS 5 years")
dev.off()


# rocCol=c("#CCCCCC","#989898","#323232")
# aucText=c()
# 
# par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
# f1_train <- cph(Surv(days_to_last_follow_up, vital_status) ~ PCBP1+PABPN1+PTPRF+DANCR+MYC, x=T, y=T, surv=T, data=rt_GSig_PFS_PANCANCER_train, time.inc = 1)
# cal1_train <- calibrate(f1_train, cmethod = 'KM', method = 'boot', u=1,m=102, B=1000)
# plot(cal1_train, type="l", xlim=c(0,1), ylim=c(0,1.05),col=rocCol[1], 
#      xlab="Predicted PFS", ylab="Observed PFS",
#      lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
# aucText=c(aucText,"1-year PFS")
# abline(0,1)
# 
# cal3_train <- calibrate(f3_train, cmethod = 'KM', method = 'boot', u=3,m=102, B=1000)
# aucText=c(aucText,"3-year PFS")
# lines(cal3_train, type="l", xlim=c(0,1), ylim=c(0,1.05),col=rocCol[2],lwd = 2)
# 
# cal5_train <- calibrate(f5_train, cmethod = 'KM', method = 'boot', u=5,m=102, B=1000)
# aucText=c(aucText,"5-year PFS")
# lines(cal5_train, type="l", xlim=c(0,1), ylim=c(0,1.05),col=rocCol[3],lwd = 2)
# 
# legend("bottomright", aucText,lwd=1,bty="n",col=rocCol)




## For testing cohort----

View(rt_GSig_PFS_PANCANCER_test)


ddlist <- datadist(rt_GSig_PFS_PANCANCER_test)
options(datadist = 'ddlist')
units(rt_GSig_PFS_PANCANCER_test$days_to_last_follow_up) <- "year"


pdf(file = "Cali_GSig_1yr_147test.pdf", width = 13, height = 8)
# 1-year calibration
f1_test <- cph(Surv(days_to_last_follow_up, vital_status) ~ PCBP1+PABPN1+PTPRF+DANCR+MYC, x=T, y=T, surv=T, data=rt_GSig_PFS_PANCANCER_test, time.inc = 1)
cal1_test <- calibrate(f1_test, cmethod = 'KM', method = 'boot', u=1,m=49, B=1000)
plot(cal1_test, xlab="Predicted PFS 1 year", ylab="Observed PFS 1 year")
dev.off()

pdf(file = "Cali_GSig_3yr_147test.pdf", width = 13, height = 8)
# 3-year calibration
f3_test <- cph(Surv(days_to_last_follow_up, vital_status) ~ PCBP1+PABPN1+PTPRF+DANCR+MYC, x=T, y=T, surv=T, data=rt_GSig_PFS_PANCANCER_test, time.inc = 3)
cal3_test <- calibrate(f3_test, cmethod = 'KM', method = 'boot', u=3,m=49, B=1000)
plot(cal3_test,xlab="Predicted PFS 3 years", ylab="Observed PFS 3 years")
dev.off()

pdf(file = "Cali_GSig_5yr_147test.pdf", width = 13, height = 8)
# 5-year calibration
f5_test <- cph(Surv(days_to_last_follow_up, vital_status) ~ PCBP1+PABPN1+PTPRF+DANCR+MYC, x=T, y=T, surv=T, data=rt_GSig_PFS_PANCANCER_test, time.inc = 5)
cal5_test <- calibrate(f5_test, cmethod = 'KM', method = 'boot', u=5,m=49, B=1000)
plot(cal5_test,xlab="Predicted PFS 5 years", ylab="Observed PFS 5 years")
dev.off()

######################################################################################################################################



##################################################################
######### Clinical info&gene signature for Cox model #############
##################################################################

## Input pre-processed clinical data----
df_clinic.genemarkers.COX <- read.table("/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Model construction/R_data_needed/df_clinic.genemarkers.COX_final.txt", header = T, sep = "\t") 
View(df_clinic.genemarkers.COX)


df_clinic.genemarkers.COX<- df_clinic.genemarkers.COX %>%
  mutate_at("sample", str_replace, "01A", "01")
df_clinic.genemarkers.COX <- df_clinic.genemarkers.COX %>%
  mutate_at("sample", str_replace, "01B", "01")
dim(df_clinic.genemarkers.COX)
##There are three duplicated samples (same tumour but A and B)
# a <- df_clinic.genemarkers.COX[df_clinic.genemarkers.COX$sample=='TCGA-HC-7740-01',]
# View(a)
df_clinic.genemarkers.COX <- df_clinic.genemarkers.COX[!duplicated(df_clinic.genemarkers.COX$sample),]
rownames(df_clinic.genemarkers.COX) <- df_clinic.genemarkers.COX$sample
View(df_clinic.genemarkers.COX)
#risk_GSig_DFS_PANCANCER.txt <- read.table("/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/DFS.PFS/risk_GSig_DFS_PANCANCER.txt", header = T, sep = "\t") 
# View(risk_GSig_DFS_PANCANCER.txt)
# rownames(risk_GSig_DFS_PANCANCER.txt) <- risk_GSig_DFS_PANCANCER.txt$id
# rownames(df_clinic.genemarkers.COX) <-df_clinic.genemarkers.COX$barcode
#NB. Columns are no longer matched, so need to re-match
#table(rownames(df_clinic.genemarkers.COX) == rownames(df_clinic.genemarkers))# check for matched or not

##
sameSample_model_train=intersect(row.names(rt_GSig_PFS_PANCANCER_train),row.names(df_clinic.genemarkers.COX)) #以'df_clinic.genemarkers' 的行名来找因为它在前而且他是最一开始的行名顺序，与‘combind_model...’一致的
sameSample_model_test=intersect(row.names(rt_GSig_PFS_PANCANCER_test),row.names(df_clinic.genemarkers.COX)) #以'df_clinic.genemarkers' 的行名来找因为它在前而且他是最一开始的行名顺序，与‘combind_model...’一致的
##

##
df_clinic.genemarkers.COX_train=df_clinic.genemarkers.COX[sameSample_model_train,]
dim(df_clinic.genemarkers.COX_train)
head(rownames(df_clinic.genemarkers.COX_train))
rt_GSig_PFS_PANCANCER_train2 <- rt_GSig_PFS_PANCANCER_train[sameSample_model_train,]
dim(rt_GSig_PFS_PANCANCER_train2)
head(rownames(rt_GSig_PFS_PANCANCER_train2))
##

##
df_clinic.genemarkers.COX_test=df_clinic.genemarkers.COX[sameSample_model_test,]
dim(df_clinic.genemarkers.COX_test)
head(rownames(df_clinic.genemarkers.COX_test))
rt_GSig_PFS_PANCANCER_test2 <- rt_GSig_PFS_PANCANCER_test[sameSample_model_test,]
dim(rt_GSig_PFS_PANCANCER_test2)
head(rownames(rt_GSig_PFS_PANCANCER_test2))
##

##Check if row names matched
table(row.names(df_clinic.genemarkers.COX_train) == row.names(rt_GSig_PFS_PANCANCER_train2))
table(row.names(df_clinic.genemarkers.COX_test) == row.names(rt_GSig_PFS_PANCANCER_test2))

## data input for uni&multi Cox regression
df_clinic.genemarkers.COX_train <- cbind(rt_GSig_PFS_PANCANCER_train2,df_clinic.genemarkers.COX_train)
View(df_clinic.genemarkers.COX_train)
dim(df_clinic.genemarkers.COX_train)
save(df_clinic.genemarkers.COX_train, file = "df_clinic.genemarkers.COX_train.rda") # Original cleaned clinical data for downstream box plot (training set)

df_clinic.genemarkers.COX_test <- cbind(rt_GSig_PFS_PANCANCER_test2,df_clinic.genemarkers.COX_test)
View(df_clinic.genemarkers.COX_test)
dim(df_clinic.genemarkers.COX_test)
save(df_clinic.genemarkers.COX_test, file = "df_clinic.genemarkers.COX_test.rda") # Original cleaned clinical data for downstream box plot (testing set)


##Clean up the df for train cohort
df_clinic.genemarkers.COX_train <- df_clinic.genemarkers.COX_train[ , -c(9:15)]#df_clinic.genemarkers.COX_train[ , -c(9:15)] #要包含之前gene signature做出来的riskScore列！
df_clinic.genemarkers.COX_train <- na.omit(df_clinic.genemarkers.COX_train)
dim(df_clinic.genemarkers.COX_train)
colnames(df_clinic.genemarkers.COX_train)
View(df_clinic.genemarkers.COX_train)
str(df_clinic.genemarkers.COX_train)
# df_clinic.genemarkers.COX_train[df_clinic.genemarkers.COX_train$vital_status == 'Alive',]$vital_status <- 0 #Alive
# df_clinic.genemarkers.COX[df_clinic.genemarkers.COX$vital_status == 'Dead',]$vital_status <- 1 #Dead
df_clinic.genemarkers.COX_train$vital_status[!df_clinic.genemarkers.COX_train$vital_status %in% c("0","1")] = NA
class(df_clinic.genemarkers.COX_train$vital_status)#class(CombindDat_model$gleason_score)
df_clinic.genemarkers.COX_train$vital_status <- as.numeric(df_clinic.genemarkers.COX_train$vital_status)#CombindDat_model$gleason_score <- as.factor(CombindDat_model$gleason_score)
table(df_clinic.genemarkers.COX_train$vital_status)

df_clinic.genemarkers.COX_train[df_clinic.genemarkers.COX_train$age_at_index == '<60',]$age_at_index <- 0 
df_clinic.genemarkers.COX_train[df_clinic.genemarkers.COX_train$age_at_index == '≥60',]$age_at_index <- 1 

df_clinic.genemarkers.COX_train[df_clinic.genemarkers.COX_train$ajcc_pathologic_t == 'T2',]$ajcc_pathologic_t <- 0 
df_clinic.genemarkers.COX_train[df_clinic.genemarkers.COX_train$ajcc_pathologic_t == 'T3-4',]$ajcc_pathologic_t <- 1 

df_clinic.genemarkers.COX_train[df_clinic.genemarkers.COX_train$ajcc_pathologic_n == 'N0',]$ajcc_pathologic_n <- 0 
df_clinic.genemarkers.COX_train[df_clinic.genemarkers.COX_train$ajcc_pathologic_n == 'N1',]$ajcc_pathologic_n <- 1 

df_clinic.genemarkers.COX_train[df_clinic.genemarkers.COX_train$gleason_score == '<8',]$gleason_score <- 0 
df_clinic.genemarkers.COX_train[df_clinic.genemarkers.COX_train$gleason_score == '≥8',]$gleason_score <- 1 

##Numerics clinic info
df_clinic.genemarkers.COX_train$age_at_index <- as.numeric(df_clinic.genemarkers.COX_train$age_at_index)
df_clinic.genemarkers.COX_train$ajcc_pathologic_t <- as.numeric(df_clinic.genemarkers.COX_train$ajcc_pathologic_t)
df_clinic.genemarkers.COX_train$ajcc_pathologic_n <- as.numeric(df_clinic.genemarkers.COX_train$ajcc_pathologic_n)
df_clinic.genemarkers.COX_train$gleason_score <- as.numeric(df_clinic.genemarkers.COX_train$gleason_score)

str(df_clinic.genemarkers.COX_train)


#Uni Cox(all variables included)----
uniTab=data.frame()
for(i in colnames(df_clinic.genemarkers.COX_train[,c(8:12)])){
  cox <- coxph(Surv(days_to_last_follow_up, vital_status) ~ df_clinic.genemarkers.COX_train[,i], data = df_clinic.genemarkers.COX_train)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(uniTab,file="uniCox_Clinc&Genesig_pfs_pan_291TRAIN.txt",sep="\t",row.names=F,quote=F)


#Multi Cox(all variables included)----
View(uniTab)
uniTab= uniTab#uniTab[as.numeric(uniTab[,"pvalue"])<0.05,]
uniTab$pvalue <- as.numeric(uniTab$pvalue)
uniTab$HR <- as.numeric(uniTab$HR)
uniTab$HR.95L <- as.numeric(uniTab$HR.95L)
uniTab$HR.95H <- as.numeric(uniTab$HR.95H)
rt1=df_clinic.genemarkers.COX_train[,c("days_to_last_follow_up","vital_status",as.vector(uniTab[,"id"]))]
View(rt1)
multiCox=coxph(Surv(days_to_last_follow_up, vital_status) ~ ., data = rt1)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file="multiCox_Clinc&Genesig_pfs_pan_291TRAIN.txt",sep="\t",row.names=F,quote=F)



##Clean up the df for test cohort
df_clinic.genemarkers.COX_test <- df_clinic.genemarkers.COX_test[ , -c(9:15)] #要包含之前gene signature做出来的riskScore列！
df_clinic.genemarkers.COX_test <- na.omit(df_clinic.genemarkers.COX_test)
dim(df_clinic.genemarkers.COX_test)
colnames(df_clinic.genemarkers.COX_test)
View(df_clinic.genemarkers.COX_test)
str(df_clinic.genemarkers.COX_test)
# df_clinic.genemarkers.COX_train[df_clinic.genemarkers.COX_train$vital_status == 'Alive',]$vital_status <- 0 #Alive
# df_clinic.genemarkers.COX[df_clinic.genemarkers.COX$vital_status == 'Dead',]$vital_status <- 1 #Dead
df_clinic.genemarkers.COX_test$vital_status[!df_clinic.genemarkers.COX_test$vital_status %in% c("0","1")] = NA
class(df_clinic.genemarkers.COX_test$vital_status)#class(CombindDat_model$gleason_score)
df_clinic.genemarkers.COX_test$vital_status <- as.numeric(df_clinic.genemarkers.COX_test$vital_status)#CombindDat_model$gleason_score <- as.factor(CombindDat_model$gleason_score)
table(df_clinic.genemarkers.COX_test$vital_status)

df_clinic.genemarkers.COX_test[df_clinic.genemarkers.COX_test$age_at_index == '<60',]$age_at_index <- 0 
df_clinic.genemarkers.COX_test[df_clinic.genemarkers.COX_test$age_at_index == '≥60',]$age_at_index <- 1 

df_clinic.genemarkers.COX_test[df_clinic.genemarkers.COX_test$ajcc_pathologic_t == 'T2',]$ajcc_pathologic_t <- 0 
df_clinic.genemarkers.COX_test[df_clinic.genemarkers.COX_test$ajcc_pathologic_t == 'T3-4',]$ajcc_pathologic_t <- 1 

df_clinic.genemarkers.COX_test[df_clinic.genemarkers.COX_test$ajcc_pathologic_n == 'N0',]$ajcc_pathologic_n <- 0 
df_clinic.genemarkers.COX_test[df_clinic.genemarkers.COX_test$ajcc_pathologic_n == 'N1',]$ajcc_pathologic_n <- 1 

df_clinic.genemarkers.COX_test[df_clinic.genemarkers.COX_test$gleason_score == '<8',]$gleason_score <- 0 
df_clinic.genemarkers.COX_test[df_clinic.genemarkers.COX_test$gleason_score == '≥8',]$gleason_score <- 1 

##Numerics clinic info
df_clinic.genemarkers.COX_test$age_at_index <- as.numeric(df_clinic.genemarkers.COX_test$age_at_index)
df_clinic.genemarkers.COX_test$ajcc_pathologic_t <- as.numeric(df_clinic.genemarkers.COX_test$ajcc_pathologic_t)
df_clinic.genemarkers.COX_test$ajcc_pathologic_n <- as.numeric(df_clinic.genemarkers.COX_test$ajcc_pathologic_n)
df_clinic.genemarkers.COX_test$gleason_score <- as.numeric(df_clinic.genemarkers.COX_test$gleason_score)

str(df_clinic.genemarkers.COX_test)




#Uni Cox(all variables included)----
uniTab=data.frame()
for(i in colnames(df_clinic.genemarkers.COX_test[,c(8:12)])){
  cox <- coxph(Surv(days_to_last_follow_up, vital_status) ~ df_clinic.genemarkers.COX_test[,i], data = df_clinic.genemarkers.COX_test)
  coxSummary = summary(cox)
  uniTab=rbind(uniTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(uniTab,file="uniCox_Clinc&Genesig_pfs_pan_125TEST.txt",sep="\t",row.names=F,quote=F)


#Multi Cox(all variables included)----
View(uniTab)
uniTab= uniTab#uniTab[as.numeric(uniTab[,"pvalue"])<0.05,]
uniTab$pvalue <- as.numeric(uniTab$pvalue)
uniTab$HR <- as.numeric(uniTab$HR)
uniTab$HR.95L <- as.numeric(uniTab$HR.95L)
uniTab$HR.95H <- as.numeric(uniTab$HR.95H)
rt2=df_clinic.genemarkers.COX_test[,c("days_to_last_follow_up","vital_status",as.vector(uniTab[,"id"]))]
View(rt2)
multiCox=coxph(Surv(days_to_last_follow_up, vital_status) ~ ., data = rt2)
multiCoxSum=summary(multiCox)
multiTab=data.frame()
multiTab=cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
multiTab=cbind(id=row.names(multiTab),multiTab)
write.table(multiTab,file="multiCox_Clinc&Genesig_pfs_pan_125TEST.txt",sep="\t",row.names=F,quote=F)

##


##################################################################
############### Forest plot for Indep prog assess ################
##################################################################

############绘制森林图函数############
bioForest=function(coxFile=null,forestFile=null,height=null,forestCol=null){
  #读取输入文件
  rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  #输出图形
  pdf(file=forestFile, width = 15,height = height)
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
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,adj=1)
  text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,adj=1,)
  
  #绘制森林图
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="black",lwd=2.5) #col="darkblue"
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol[1], forestCol[2])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
  axis(1)
  dev.off()
}

bioForest(coxFile="uniCox_Clinc&Genesig_pfs_pan_291TRAIN.txt",forestFile="uniForest_train.pdf",height=3.5,forestCol=c("black","black"))
bioForest(coxFile="multiCox_Clinc&Genesig_pfs_pan_291TRAIN.txt",forestFile="multiForest_train.pdf",height=3.5,forestCol=c("black","black"))
bioForest(coxFile="uniCox_Clinc&Genesig_pfs_pan_125TEST.txt",forestFile="uniForest_test.pdf",height=3.5,forestCol=c("black","black"))
bioForest(coxFile="multiCox_Clinc&Genesig_pfs_pan_125TEST.txt",forestFile="multiForest_test.pdf",height=3.5,forestCol=c("black","black"))




##################################################################
######## Nomogram (riskscore+significant clinical info) ##########
##################################################################
library(rms)
library(Hmisc)
library(lattice)
library(Formula)
library(ggplot2)
library(foreign)

#lnomo_train <- load('/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/DFS.PFS/PFS_edgeR/split_kym/final_v2/nomo_train.RData')
nomo_train <- df_clinic.genemarkers.COX_train # Filtered one
nomo_train <- nomo_train[ , -c(3:7,9:11)]
View(nomo_train)
colnames(nomo_train) <- c('futime','fustat','Five_gene_riskscore','Gleason_score')
str(nomo_train)
#nomo_train$Gleason_score <- as.factor(nomo_train$Gleason_score) #factor or numerics doesnt matter

#data pack
dd <- datadist(nomo_train)
options(datadist="dd")
summary(nomo_train$futime)

#mul cox model
f <- cph(Surv(futime, fustat) ~ Five_gene_riskscore+Gleason_score, x=T, y=T, surv=T, data=nomo_train, time.inc=1) #对于nomogram好像加不加‘time.inc=1’都一样; sxzxw 加了
print(f)
surv <- Survival(f)



#nomogram set up
nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(3, x), function(x) surv(5, x)), 
                lp=F, funlabel=c("1-year PFS", "3-year PFS", "5-year PFS"), 
                maxscale=100, 
                #fun.at=c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1,0.05)
                fun.at=c(0.95, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3,0.2,0.1,0.05)# paper use 0.9 not 0.99 a lot
                )  

print(nom)
#install.packages("nomogramEx")
library(nomogramEx)
nomogramEx(nomo=nom,np=2,digit = 9)

# #nomogram visualisation opt1
# pdf(file="nomogram_291train.pdf",height=6,width=10)
# plot(nom)
# dev.off()


#nomogram visualisation opt2
pdf(file="nomogram_291train_2.0.pdf",height=6,width=10)
plot(nom, xfrac=.35,varname.label=TRUE, varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE,cex.var = 1.6,cex.axis = 1.05,lwd=5,
     label.every = 1,col.grid = gray(c(0.8, 0.95)) )
dev.off()





## Calibration plot For training cohort----
View(df_clinic.genemarkers.COX_train)

#data packing
ddlist <- datadist(df_clinic.genemarkers.COX_train)
options(datadist = 'ddlist')
units(df_clinic.genemarkers.COX_train$days_to_last_follow_up) <- "year"

pdf(file = "Cali_GSigPlusGleason_1yr_291train.pdf", width = 13, height = 8)
# 1-year calibration
f1_train <- cph(Surv(days_to_last_follow_up, vital_status) ~ riskScore_Gsig_PFS_PANCANCER_train+gleason_score, x=T, y=T, surv=T, data=df_clinic.genemarkers.COX_train, time.inc = 1)
cal1_train <- calibrate(f1_train, cmethod = 'KM', method = 'boot', u=1,m=97, B=1000)
plot(cal1_train, xlab="Nomogram-predicted PFS 1 year", ylab="Observed PFS 1 year")
dev.off()

pdf(file = "Cali_GSigPlusGleason_3yr_291train.pdf", width = 13, height = 8)
# 3-year calibration
f3_train <- cph(Surv(days_to_last_follow_up, vital_status) ~ riskScore_Gsig_PFS_PANCANCER_train+gleason_score, x=T, y=T, surv=T, data=df_clinic.genemarkers.COX_train, time.inc = 3)
cal3_train <- calibrate(f3_train, cmethod = 'KM', method = 'boot', u=3,m=97, B=1000)
plot(cal3_train,xlab="Nomogram-predicted PFS 3 years", ylab="Observed PFS 3 years")
dev.off()

pdf(file = "Cali_GSigPlusGleason_5yr_291train.pdf", width = 13, height = 8)
# 5-year calibration
f5_train <- cph(Surv(days_to_last_follow_up, vital_status) ~ riskScore_Gsig_PFS_PANCANCER_train+gleason_score, x=T, y=T, surv=T, data=df_clinic.genemarkers.COX_train, time.inc = 5)
cal5_train <- calibrate(f5_train, cmethod = 'KM', method = 'boot', u=5,m=97, B=1000)
plot(cal5_train,xlab="Nomogram-predicted PFS 5 years", ylab="Observed PFS 5 years")
dev.off()



# Combined version!!!
pdf("calibration_135compare_291train.pdf",width = 8,height = 8)
plot(cal1_train,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced PFS probability",ylab = "Observed PFS probability",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1_train[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal3_train,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal3_train[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

plot(cal5_train,lwd = 2,lty = 0,errbar.col = c("#2A960C"),
     xlim = c(0,1),ylim= c(0,1),col = c("#2A960C"),add = T)
lines(cal5_train[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2A960C"), pch = 16)


abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","3-year","5-year"), 
       col =c("#2166AC","#B2182B","#2A960C"), 
       lwd = 2,
       cex = 1.2,
       bty = "n")
dev.off()



## Calibration plot For testing cohort----
View(df_clinic.genemarkers.COX_test)

#data packing
ddlist <- datadist(df_clinic.genemarkers.COX_test)
options(datadist = 'ddlist')
units(df_clinic.genemarkers.COX_test$days_to_last_follow_up) <- "year"

pdf(file = "Cali_GSigPlusGleason_1yr_125test.pdf", width = 13, height = 8)
# 1-year calibration
f1_test <- cph(Surv(days_to_last_follow_up, vital_status) ~ riskScore_Gsig_PFS_PANCANCER_test+gleason_score, x=T, y=T, surv=T, data=df_clinic.genemarkers.COX_test, time.inc = 1)
cal1_test <- calibrate(f1_test, cmethod = 'KM', method = 'boot', u=1,m=41, B=1000)
plot(cal1_test, xlab="Nomogram-predicted PFS 1 year", ylab="Observed PFS 1 year")
dev.off()

pdf(file = "Cali_GSigPlusGleason_3yr_125test.pdf", width = 13, height = 8)
# 3-year calibration
f3_test <- cph(Surv(days_to_last_follow_up, vital_status) ~ riskScore_Gsig_PFS_PANCANCER_test+gleason_score, x=T, y=T, surv=T, data=df_clinic.genemarkers.COX_test, time.inc = 3)
cal3_test <- calibrate(f3_test, cmethod = 'KM', method = 'boot', u=3,m=41, B=1000)
plot(cal3_test,xlab="Nomogram-predicted PFS 3 years", ylab="Observed PFS 3 years")
dev.off()

pdf(file = "Cali_GSigPlusGleason_5yr_125test.pdf", width = 13, height = 8)
# 5-year calibration
f5_test <- cph(Surv(days_to_last_follow_up, vital_status) ~ riskScore_Gsig_PFS_PANCANCER_test+gleason_score, x=T, y=T, surv=T, data=df_clinic.genemarkers.COX_test, time.inc = 5)
cal5_test <- calibrate(f5_test, cmethod = 'KM', method = 'boot', u=5,m=41, B=1000)
plot(cal5_test,xlab="Nomogram-predicted PFS 5 years", ylab="Observed PFS 5 years")
dev.off()

# Combined version!!!
pdf("calibration_135compare_125test.pdf",width = 8,height = 8)
plot(cal1_test,lwd = 2,lty = 0,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced PFS probability",ylab = "Observed PFS probability",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1_test[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal3_test,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal3_test[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)

plot(cal5_test,lwd = 2,lty = 0,errbar.col = c("#2A960C"),
     xlim = c(0,1),ylim= c(0,1),col = c("#2A960C"),add = T)
lines(cal5_test[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2A960C"), pch = 16)


abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","3-year","5-year"), 
       col =c("#2166AC","#B2182B","#2A960C"), 
       lwd = 2,
       cex = 1.2,
       bty = "n")
dev.off()



##################################################################
########### Multivariate COX for (GS and Gleason score) #########
##################################################################
# Multivariate COX for (GSig riskscore and Gleason score)----
#!!!Questions: Need to run the ORIGINAL sample data or the clinical NA filtered one for construction of gleason score related model?
df_clinic.genemarkers.COX_train$riskScore <- df_clinic.genemarkers.COX_train$riskScore_Gsig_PFS_PANCANCER_train
View(df_clinic.genemarkers.COX_train)
df_clinic.genemarkers.COX_test$riskScore <- df_clinic.genemarkers.COX_test$riskScore_Gsig_PFS_PANCANCER_test
View(df_clinic.genemarkers.COX_test)

#Model construction
cox_Mul2 <- coxph(Surv(days_to_last_follow_up,vital_status) ~ gleason_score+riskScore , data = df_clinic.genemarkers.COX_train) #age_at_index+ajcc_pathologic_t+ajcc_pathologic_n+gleason_score+MRPL24+MYC+HNRNPL+ELAVL1+IMPDH2#coxph(Surv(days_to_last_follow_up,vital_status) ~ ., data = CombindDat_model)
cox_Mul2
#cox_Mul2 <- step(cox_Mul2, direction = "both") # 几乎什么情况都是gleason 被选出; 上边一步中直接写出这步得到变量也可以的，结果一样
riskScore_GSigPlusGS_train <- predict(cox_Mul2, type = "risk", newdata = df_clinic.genemarkers.COX_train)
riskScore_GSigPlusGS_test <- predict(cox_Mul2, type = "risk", newdata = df_clinic.genemarkers.COX_test)

summary <- summary(cox_Mul2)

library(Hmisc)
1-with(df_clinic.genemarkers.COX_train, rcorr.cens(riskScore_GSigPlusGS_train, Surv(days_to_last_follow_up,vital_status)))[1]
1-with(df_clinic.genemarkers.COX_test, rcorr.cens(riskScore_GSigPlusGS_test, Surv(days_to_last_follow_up,vital_status)))[1]

outCol <- c("days_to_last_follow_up","vital_status","gleason_score","riskScore")
#outCol <- c("days_to_last_follow_up","vital_status",coxGene)
View(outCol)
##Now, we use optiml cut-off to classify riskscore on the Train and Test set now! no more on the same training set as before
df_clinic.genemarkers.COX_train$risk <- riskScore_GSigPlusGS_train
View(df_clinic.genemarkers.COX_train)
df_clinic.genemarkers.COX_test$risk <- riskScore_GSigPlusGS_test
View(df_clinic.genemarkers.COX_test)




###################### Optimal cut-ff approach for survival analysis of calculated riskScore (TRAINING cohort) #######################
surv.cut_PRAD <- surv_cutpoint(
  df_clinic.genemarkers.COX_train,
  time = "days_to_last_follow_up",
  event = "vital_status",
  variables = "risk"
)

summary(surv.cut_PRAD)
surv.cat_PRAD <- surv_categorize(surv.cut_PRAD)
View(surv.cat_PRAD)
diff <- survdiff(Surv(days_to_last_follow_up,vital_status) ~ risk, data = surv.cat_PRAD)
pValue <- 1-pchisq(diff$chisq, df=1)
pValue <- signif(pValue, 4)
pValue <- format(pValue, scientific=FALSE)

fit <- survfit(Surv(days_to_last_follow_up,vital_status) ~ risk, data = surv.cat_PRAD)
summary(fit) #check for 5-year survival rate

###################### Plot: Overall/progression-free survival for GSig ########################
# tiff(file = "survival_MulCox_KMplot_MYC&ELAVL1_fpkm.tiff",
#      width = 14,
#      height = 14,
#      units = "cm",
#      compression = "lzw",
#      bg = "white",
#      res = 600)
#
#PDF
# pdf(file = "survival_MulCox_KMplot_GSigPlusGleason_edgeR_pfs_pan_256train.pdf", width = 8, height = 8)
# 
# 
# 
plot(fit,
     lwd = 2,
     col = c("#323232","#CCCCCC"),#col = c("red","blue"),
     xlab = "Time (year)",
     mark.time = T,
     ylab = "Progression-free survival"
     # main = paste("Progression-free survival curve (p=", pValue ,")", sep = "")
)
legend("topright",
       c("high risk", "low risk"),
       lwd = 2,
       col = c("#323232","#CCCCCC"))
legend("bottomright",
       paste("Logrank p=",pValue, sep = ""))
# 
# dev.off()
####



#Survival curve for MulCox Option2 (with risk table)----
#KM analysis and survival curve based on the high/low risk
fit_PRAD <- survfit(Surv(days_to_last_follow_up,vital_status) ~ risk,
                    data = surv.cat_PRAD)

pdf(file = "survival_MulCox_KMplot_GSigPlusGleason_edgeR_pfs_pan_291train.pdf", width = 8, height = 8)

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


#riskScore table export for training set----
## re-define the risk to make sure to include corretly in the final 'risk' table
risk <- surv.cat_PRAD$risk
#risk <- as.vector(ifelse(riskScore>median(riskScore), "high","low")) #No more MEDIAN！！！
View(df_clinic.genemarkers.COX_train)
##write out 'risk' table
write.table(cbind(id = rownames(cbind(df_clinic.genemarkers.COX_train[ , outCol], riskScore_GSigPlusGS_train, risk)), cbind(df_clinic.genemarkers.COX_train[ , outCol],riskScore_GSigPlusGS_train,risk)),
            file = "risk_GSigPlusGleason_PFS_PANCANCER_291train.txt",
            sep = "\t",
            quote = F,
            row.names = F)

rt_GSigPlusGleason_PFS_PANCANCER_train <- read.table("risk_GSigPlusGleason_PFS_PANCANCER_291train.txt", header = T, sep = "\t", check.names = F, row.names = 1)
View(rt_GSigPlusGleason_PFS_PANCANCER_train)

######################################################################################################################################



###################### Optimal cut-ff approach for survival analysis of calculated riskScore (TESTING cohort) #######################

surv.cut_PRAD <- surv_cutpoint(
  df_clinic.genemarkers.COX_test,
  time = "days_to_last_follow_up",
  event = "vital_status",
  variables = "risk"
)

summary(surv.cut_PRAD)
surv.cat_PRAD <- surv_categorize(surv.cut_PRAD)
View(surv.cat_PRAD)
diff <- survdiff(Surv(days_to_last_follow_up,vital_status) ~ risk, data = surv.cat_PRAD)
pValue <- 1-pchisq(diff$chisq, df=1)
pValue <- signif(pValue, 4)
pValue <- format(pValue, scientific=FALSE)

fit <- survfit(Surv(days_to_last_follow_up,vital_status) ~ risk, data = surv.cat_PRAD)
summary(fit) #check for 5-year survival rate

###################### Plot: Overall/progression-free survival for GSig ########################
# tiff(file = "survival_MulCox_KMplot_MYC&ELAVL1_fpkm.tiff",
#      width = 14,
#      height = 14,
#      units = "cm",
#      compression = "lzw",
#      bg = "white",
#      res = 600)
#
#PDF
# pdf(file = "survival_MulCox_KMplot_GSigPlusGleason_edgeR_pfs_pan_168test.pdf", width = 8, height = 8)
# 
# 
# 
plot(fit,
     lwd = 2,
     col = c("#323232","#CCCCCC"),#col = c("red","blue"),
     xlab = "Time (year)",
     mark.time = T,
     ylab = "Progression-free survival"
     # main = paste("Progression-free survival curve (p=", pValue ,")", sep = "")
)
legend("topright",
       c("high risk", "low risk"),
       lwd = 2,
       col = c("#323232","#CCCCCC"))
legend("bottomright",
       paste("Logrank p=",pValue, sep = ""))
# 
# dev.off()
####


#Survival curve for MulCox Option2 (with risk table)----
#KM analysis and survival curve based on the high/low risk
fit_PRAD <- survfit(Surv(days_to_last_follow_up,vital_status) ~ risk,
                    data = surv.cat_PRAD)

pdf(file = "survival_MulCox_KMplot_GSigPlusGleason_edgeR_pfs_pan_125test.pdf", width = 8, height = 8)

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



#riskScore table export for testing set----
## re-define the risk to make sure to include corretly in the final 'risk' table
risk <- surv.cat_PRAD$risk
#risk <- as.vector(ifelse(riskScore>median(riskScore), "high","low")) #No more MEDIAN！！！
View(df_clinic.genemarkers.COX_test)
##write out 'risk' table
write.table(cbind(id = rownames(cbind(df_clinic.genemarkers.COX_test[ , outCol], riskScore_GSigPlusGS_test, risk)), cbind(df_clinic.genemarkers.COX_test[ , outCol],riskScore_GSigPlusGS_test,risk)),
            file = "risk_GSigPlusGleason_PFS_PANCANCER_125test.txt",
            sep = "\t",
            quote = F,
            row.names = F)

rt_GSigPlusGleason_PFS_PANCANCER_test <- read.table("risk_GSigPlusGleason_PFS_PANCANCER_125test.txt", header = T, sep = "\t", check.names = F, row.names = 1)
View(rt_GSigPlusGleason_PFS_PANCANCER_test)

######################################################################################################################################


###################################################### Time-dependent ROC analysis ###################################################
##From: https://www.jianshu.com/p/b5e96158059b?utm_campaign=hugo&utm_medium=reader_share&utm_content=note

## For training set----
View(rt_GSigPlusGleason_PFS_PANCANCER_train)
colnames(rt_GSigPlusGleason_PFS_PANCANCER_train)
result <-with(rt_GSigPlusGleason_PFS_PANCANCER_train, timeROC( T=days_to_last_follow_up,
                                                               delta=vital_status,
                                                               marker=riskScore_GSigPlusGS_train,#marker=riskScore_Gsig_PFS_PANCANCER_train_modfy,
                                                               cause=1,
                                                               times=c(1,3,5),
                                                               iid = TRUE))
identical(c(result$TP[,1],result$TP[,2],result$TP[,3]),as.numeric(result$TP))
dat = data.frame(fpr = as.numeric(result$FP),
                 tpr = as.numeric(result$TP),
                 time = rep(as.factor(c(1,3,5)),each = nrow(result$TP)))



# #TIFF
# tiff(file = "MulCox_ROC_allyrs_2Genesig_fpkm.tiff",
#      width = 20,
#      height = 20,
#      units = "cm",
#      compression = "lzw",
#      bg = "white",
#      res = 600)

#PDF for training set
pdf(file = "ROC_135yrs_GenesigPlusGleason_291TCGAvalid.pdf", width = 8, height = 8)


par(oma = c(0.5,1,0,1), font.lab = 1.5, font.axis = 1.5)

ggplot() +
  geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) +
  scale_color_manual(name = NULL,values = c("#323232", "#989898", "#CCCCCC"),
                     labels = paste0("AUC of ",c(1,3,5),"-year survival: ",
                                     format(round(result$AUC,3),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()+
  ggtitle("ROC curve of gene signature and Gleason score for PFS in TCGA-PRAD training cohort")+
  theme(plot.title = element_text(hjust = 0.5,colour = "black",face = "bold.italic",size = 10))

dev.off()

## Alt ROC2
# library(survivalROC)
# 
# 
# pdf(file="ROC.pdf",width=6,height=6)
# par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
# roc=survivalROC(Stime=rt_GSig_PFS_PANCANCER_train$days_to_last_follow_up, status=rt_GSig_PFS_PANCANCER_train$vital_status, marker = rt_GSig_PFS_PANCANCER_train$riskScore_Gsig_PFS_PANCANCER_train, 
#                 predict.time =5, method="KM")
# plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
#      xlab="False positive rate", ylab="True positive rate",
#      main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
#      lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
# abline(0,1)
# dev.off()




## For testing set----
View(rt_GSigPlusGleason_PFS_PANCANCER_test)
colnames(rt_GSigPlusGleason_PFS_PANCANCER_test)
result <-with(rt_GSigPlusGleason_PFS_PANCANCER_test, timeROC( T=days_to_last_follow_up,
                                                              delta=vital_status,
                                                              marker=riskScore_GSigPlusGS_test,#marker=riskScore_Gsig_PFS_PANCANCER_test_modfy,
                                                              cause=1,
                                                              times=c(1,3,5),
                                                              iid = TRUE))
identical(c(result$TP[,1],result$TP[,2],result$TP[,3]),as.numeric(result$TP))
dat = data.frame(fpr = as.numeric(result$FP),
                 tpr = as.numeric(result$TP),
                 time = rep(as.factor(c(1,3,5)),each = nrow(result$TP)))

#PDF for testing set
pdf(file = "ROC_135yrs_GenesigPlusGleason_125TCGAvalid.pdf", width = 8, height = 8)


par(oma = c(0.5,1,0,1), font.lab = 1.5, font.axis = 1.5)

ggplot() +
  geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) +
  scale_color_manual(name = NULL,values = c("#323232", "#989898", "#CCCCCC"),
                     labels = paste0("AUC of ",c(1,3,5),"-year survival: ",
                                     format(round(result$AUC,3),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()+
  ggtitle("ROC curve of gene signature and Gleason score for PFS in TCGA-PRAD testing cohort")+
  theme(plot.title = element_text(hjust = 0.5,colour = "black",face = "bold.italic",size = 10))

dev.off()
##

# Alt ROC2
# library(survivalROC)
# 
# 
# pdf(file="ROC.pdf",width=6,height=6)
# par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
# roc=survivalROC(Stime=rt_GSig_PFS_PANCANCER_test$days_to_last_follow_up, status=rt_GSig_PFS_PANCANCER_test$vital_status, marker = rt_GSig_PFS_PANCANCER_test$riskScore_Gsig_PFS_PANCANCER_test,
#                 predict.time =3, method="KM") #predict 3 year
# plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red',
#      xlab="False positive rate", ylab="True positive rate",
#      main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
#      lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
# abline(0,1)
# dev.off()


#c-index， similar to AUC， often used in COX evaluating cancer patients prognostic model's accuracy----
#BiocManager::install("survcomp")
library(survcomp)
# For training set
cindex <- concordance.index(x = rt_GSigPlusGleason_PFS_PANCANCER_train$riskScore_GSigPlusGS_train,
                            surv.time = rt_GSigPlusGleason_PFS_PANCANCER_train$days_to_last_follow_up,
                            surv.event = rt_GSigPlusGleason_PFS_PANCANCER_train$vital_status,
                            method = "noether")

cindex$c.index
cindex$se
cindex$lower
cindex$upper
cindex$p.value

#For testing set
cindex <- concordance.index(x = rt_GSigPlusGleason_PFS_PANCANCER_test$riskScore_GSigPlusGS_test,
                            surv.time = rt_GSigPlusGleason_PFS_PANCANCER_test$days_to_last_follow_up,
                            surv.event = rt_GSigPlusGleason_PFS_PANCANCER_test$vital_status,
                            method = "noether")

cindex$c.index
cindex$se
cindex$lower
cindex$upper
cindex$p.value

######################################################################################################################################




##################################################################
########### Multivariate COX for (Gleason score ONLY) ############
##################################################################
# Multivariate COX for (GSig riskscore and Gleason score)----

dim(df_clinic.genemarkers.COX_train)

#Model construction
cox_Mul3 <- coxph(Surv(days_to_last_follow_up,vital_status) ~ gleason_score, data = df_clinic.genemarkers.COX_train) #age_at_index+ajcc_pathologic_t+ajcc_pathologic_n+gleason_score+MRPL24+MYC+HNRNPL+ELAVL1+IMPDH2#coxph(Surv(days_to_last_follow_up,vital_status) ~ ., data = CombindDat_model)
#cox_Mul2 <- coxph(Surv(days_to_last_follow_up,vital_status) ~ gleason_score+PABPN1+TRPM4+PCBP1 , data = df_clinic.genemarkers.COX_train) #age_at_index+ajcc_pathologic_t+ajcc_pathologic_n+gleason_score+MRPL24+MYC+HNRNPL+ELAVL1+IMPDH2#coxph(Surv(days_to_last_follow_up,vital_status) ~ ., data = CombindDat_model)



#cox_Mul2 <- step(cox_Mul2, direction = "both") # 几乎什么情况都是gleason 被选出; 上边一步中直接写出这步得到变量也可以的，结果一样
riskScore_GS_train <- predict(cox_Mul3, type = "risk", newdata = df_clinic.genemarkers.COX_train)
riskScore_GS_test <- predict(cox_Mul3, type = "risk", newdata = df_clinic.genemarkers.COX_test)

summary <- summary(cox_Mul3)
# coxGene <- rownames(summary$coefficients)
# coxGene <- gsub("`","",coxGene)
library(Hmisc)
1-with(df_clinic.genemarkers.COX_train, rcorr.cens(riskScore_GS_train, Surv(days_to_last_follow_up,vital_status)))[1]
1-with(df_clinic.genemarkers.COX_test, rcorr.cens(riskScore_GS_test, Surv(days_to_last_follow_up,vital_status)))[1]



# coxGene <- rownames(summary$coefficients)
# coxGene <- gsub("`","",coxGene)
#outCol <- c("days_to_last_follow_up","vital_status",coxGene)
outCol <- c("days_to_last_follow_up","vital_status","gleason_score")
View(outCol)
##Now, we use optiml cut-off to classify riskscore on the Train and Test set now! no more on the same training set as before
df_clinic.genemarkers.COX_train$risk <- riskScore_GS_train
View(df_clinic.genemarkers.COX_train)
df_clinic.genemarkers.COX_test$risk <- riskScore_GS_test
View(df_clinic.genemarkers.COX_test)




###################### Optimal cut-ff approach for survival analysis of calculated riskScore (TRAINING cohort) #######################
surv.cut_PRAD <- surv_cutpoint(
  df_clinic.genemarkers.COX_train,
  time = "days_to_last_follow_up",
  event = "vital_status",
  variables = "risk"
)

summary(surv.cut_PRAD)
surv.cat_PRAD <- surv_categorize(surv.cut_PRAD)
View(surv.cat_PRAD)
diff <- survdiff(Surv(days_to_last_follow_up,vital_status) ~ risk, data = surv.cat_PRAD)
pValue <- 1-pchisq(diff$chisq, df=1)
pValue <- signif(pValue, 4)
pValue <- format(pValue, scientific=FALSE)

fit <- survfit(Surv(days_to_last_follow_up,vital_status) ~ risk, data = surv.cat_PRAD)
summary(fit) #check for 5-year survival rate

###################### Plot: Overall/progression-free survival for GSig ########################
# tiff(file = "survival_MulCox_KMplot_MYC&ELAVL1_fpkm.tiff",
#      width = 14,
#      height = 14,
#      units = "cm",
#      compression = "lzw",
#      bg = "white",
#      res = 600)
#
#PDF
# pdf(file = "survival_MulCox_KMplot_Gleason_edgeR_pfs_pan_222train.pdf", width = 8, height = 8)
# 
# 
# 
# plot(fit,
#      lwd = 2,
#      col = c("#323232","#CCCCCC"),#col = c("red","blue"),
#      xlab = "Time (year)",
#      mark.time = T,
#      ylab = "Progression-free survival"
#      # main = paste("Progression-free survival curve (p=", pValue ,")", sep = "")
# )
# legend("topright",
#        c("high risk", "low risk"),
#        lwd = 2,
#        col = c("#323232","#CCCCCC"))
# legend("bottomright",
#        paste("Logrank p=",pValue, sep = ""))
# 
# dev.off()
####



#riskScore table export for training set----
## re-define the risk to make sure to include corretly in the final 'risk' table
risk <- surv.cat_PRAD$risk
#risk <- as.vector(ifelse(riskScore>median(riskScore), "high","low")) #No more MEDIAN！！！
View(df_clinic.genemarkers.COX_train)
##write out 'risk' table
write.table(cbind(id = rownames(cbind(df_clinic.genemarkers.COX_train[ , outCol], riskScore_GS_train, risk)), cbind(df_clinic.genemarkers.COX_train[ , outCol],riskScore_GS_train,risk)),
            file = "risk_Gleason_PFS_PANCANCER_291train.txt",
            sep = "\t",
            quote = F,
            row.names = F)


rt_Gleason_PFS_PANCANCER_train <- read.table("risk_Gleason_PFS_PANCANCER_291train.txt", header = T, sep = "\t", check.names = F, row.names = 1)
View(rt_Gleason_PFS_PANCANCER_train)

######################################################################################################################################



###################### Optimal cut-ff approach for survival analysis of calculated riskScore (TESTING cohort) #######################

surv.cut_PRAD <- surv_cutpoint(
  df_clinic.genemarkers.COX_test,
  time = "days_to_last_follow_up",
  event = "vital_status",
  variables = "risk"
)

summary(surv.cut_PRAD)
surv.cat_PRAD <- surv_categorize(surv.cut_PRAD)
View(surv.cat_PRAD)
diff <- survdiff(Surv(days_to_last_follow_up,vital_status) ~ risk, data = surv.cat_PRAD)
pValue <- 1-pchisq(diff$chisq, df=1)
pValue <- signif(pValue, 4)
pValue <- format(pValue, scientific=FALSE)

fit <- survfit(Surv(days_to_last_follow_up,vital_status) ~ risk, data = surv.cat_PRAD)
summary(fit) #check for 5-year survival rate

###################### Plot: Overall/progression-free survival for GSig ########################
# tiff(file = "survival_MulCox_KMplot_MYC&ELAVL1_fpkm.tiff",
#      width = 14,
#      height = 14,
#      units = "cm",
#      compression = "lzw",
#      bg = "white",
#      res = 600)
#
#PDF
# pdf(file = "survival_MulCox_KMplot_GSigPlusGleason_edgeR_pfs_pan_168test.pdf", width = 8, height = 8)
# 
# 
# 
# plot(fit,
#      lwd = 2,
#      col = c("#323232","#CCCCCC"),#col = c("red","blue"),
#      xlab = "Time (year)",
#      mark.time = T,
#      ylab = "Progression-free survival"
#      # main = paste("Progression-free survival curve (p=", pValue ,")", sep = "")
# )
# legend("topright",
#        c("high risk", "low risk"),
#        lwd = 2,
#        col = c("#323232","#CCCCCC"))
# legend("bottomright",
#        paste("Logrank p=",pValue, sep = ""))
# 
# dev.off()
####




#riskScore table export for testing set----
## re-define the risk to make sure to include corretly in the final 'risk' table
risk <- surv.cat_PRAD$risk
#risk <- as.vector(ifelse(riskScore>median(riskScore), "high","low")) #No more MEDIAN！！！
View(df_clinic.genemarkers.COX_test)
##write out 'risk' table
write.table(cbind(id = rownames(cbind(df_clinic.genemarkers.COX_test[ , outCol], riskScore_GS_test, risk)), cbind(df_clinic.genemarkers.COX_test[ , outCol],riskScore_GS_test,risk)),
            file = "risk_Gleason_PFS_PANCANCER_125test.txt",
            sep = "\t",
            quote = F,
            row.names = F)

rt_Gleason_PFS_PANCANCER_test <- read.table("risk_Gleason_PFS_PANCANCER_125test.txt", header = T, sep = "\t", check.names = F, row.names = 1)
View(rt_Gleason_PFS_PANCANCER_test)

######################################################################################################################################


###################################################### Time-dependent ROC analysis ###################################################
##From: https://www.jianshu.com/p/b5e96158059b?utm_campaign=hugo&utm_medium=reader_share&utm_content=note

## For training set----
View(rt_Gleason_PFS_PANCANCER_train)
colnames(rt_Gleason_PFS_PANCANCER_train)
result <-with(rt_Gleason_PFS_PANCANCER_train, timeROC( T=days_to_last_follow_up,
                                                       delta=vital_status,
                                                       marker=riskScore_GS_train,#marker=riskScore_Gsig_PFS_PANCANCER_train_modfy,
                                                       cause=1,
                                                       times=c(1,3,5),
                                                       iid = TRUE))
identical(c(result$TP[,1],result$TP[,2],result$TP[,3]),as.numeric(result$TP))
dat = data.frame(fpr = as.numeric(result$FP),
                 tpr = as.numeric(result$TP),
                 time = rep(as.factor(c(1,3,5)),each = nrow(result$TP)))



# #TIFF
# tiff(file = "MulCox_ROC_allyrs_2Genesig_fpkm.tiff",
#      width = 20,
#      height = 20,
#      units = "cm",
#      compression = "lzw",
#      bg = "white",
#      res = 600)

#PDF for training set
pdf(file = "ROC_135yrs_Gleason_291TCGAvalid.pdf", width = 8, height = 8)


par(oma = c(0.5,1,0,1), font.lab = 1.5, font.axis = 1.5)

ggplot() +
  geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) +
  scale_color_manual(name = NULL,values = c("#323232", "#989898", "#CCCCCC"),
                     labels = paste0("AUC of ",c(1,3,5),"-year survival: ",
                                     format(round(result$AUC,3),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()+
  ggtitle("ROC curve of Gleason score for PFS in TCGA-PRAD training cohort")+
  theme(plot.title = element_text(hjust = 0.5,colour = "black",face = "bold.italic",size = 10))

dev.off()

## Alt ROC2
# library(survivalROC)
# 
# 
# pdf(file="ROC.pdf",width=6,height=6)
# par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
# roc=survivalROC(Stime=rt_GSig_PFS_PANCANCER_train$days_to_last_follow_up, status=rt_GSig_PFS_PANCANCER_train$vital_status, marker = rt_GSig_PFS_PANCANCER_train$riskScore_Gsig_PFS_PANCANCER_train, 
#                 predict.time =5, method="KM")
# plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
#      xlab="False positive rate", ylab="True positive rate",
#      main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
#      lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
# abline(0,1)
# dev.off()




## For testing set----
View(rt_Gleason_PFS_PANCANCER_test)
colnames(rt_Gleason_PFS_PANCANCER_test)
result <-with(rt_Gleason_PFS_PANCANCER_test, timeROC( T=days_to_last_follow_up,
                                                      delta=vital_status,
                                                      marker=riskScore_GS_test,#marker=riskScore_Gsig_PFS_PANCANCER_test_modfy,
                                                      cause=1,
                                                      times=c(1,3,5),
                                                      iid = TRUE))
identical(c(result$TP[,1],result$TP[,2],result$TP[,3]),as.numeric(result$TP))
dat = data.frame(fpr = as.numeric(result$FP),
                 tpr = as.numeric(result$TP),
                 time = rep(as.factor(c(1,3,5)),each = nrow(result$TP)))

#PDF for testing set
pdf(file = "ROC_135yrs_Gleason_125TCGAvalid.pdf", width = 8, height = 8)


par(oma = c(0.5,1,0,1), font.lab = 1.5, font.axis = 1.5)

ggplot() +
  geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) +
  scale_color_manual(name = NULL,values = c("#323232", "#989898", "#CCCCCC"),
                     labels = paste0("AUC of ",c(1,3,5),"-year survival: ",
                                     format(round(result$AUC,3),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()+
  ggtitle("ROC curve of Gleason score for PFS in TCGA-PRAD testing cohort")+
  theme(plot.title = element_text(hjust = 0.5,colour = "black",face = "bold.italic",size = 10))

dev.off()
##

# Alt ROC2
# library(survivalROC)
# 
# 
# pdf(file="ROC.pdf",width=6,height=6)
# par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
# roc=survivalROC(Stime=rt_GSig_PFS_PANCANCER_test$days_to_last_follow_up, status=rt_GSig_PFS_PANCANCER_test$vital_status, marker = rt_GSig_PFS_PANCANCER_test$riskScore_Gsig_PFS_PANCANCER_test,
#                 predict.time =3, method="KM") #predict 3 year
# plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red',
#      xlab="False positive rate", ylab="True positive rate",
#      main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
#      lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
# abline(0,1)
# dev.off()


#c-index， similar to AUC， often used in COX evaluating cancer patients prognostic model's accuracy----
#BiocManager::install("survcomp")
library(survcomp)

# For training set
cindex <- concordance.index(x = rt_Gleason_PFS_PANCANCER_train$riskScore_GS_train,
                            surv.time = rt_Gleason_PFS_PANCANCER_train$days_to_last_follow_up,
                            surv.event = rt_Gleason_PFS_PANCANCER_train$vital_status,
                            method = "noether")

cindex$c.index
cindex$se
cindex$lower
cindex$upper
cindex$p.value

# For testing set
cindex <- concordance.index(x = rt_Gleason_PFS_PANCANCER_test$riskScore_GS_test,
                            surv.time = rt_Gleason_PFS_PANCANCER_test$days_to_last_follow_up,
                            surv.event = rt_Gleason_PFS_PANCANCER_test$vital_status,
                            method = "noether")

cindex$c.index
cindex$se
cindex$lower
cindex$upper
cindex$p.value

######################################################################################################################################



#####################################################################################
#### 1/3/5-year ROC curves for all three types models (B&W) using 'survivalROC'######
#####################################################################################
#NEED to CHANGE YEAR MANNUALLY
## For training cohort----

library(survivalROC)

## run together！！！
rocCol=c("#CCCCCC","#989898","#323232")
aucText=c()

#ROC curve
##TIFF
# tiff(file="ROCs_all_3yr.tiff",
#      width=18,
#      height=18,
#      units ="cm",
#      compression="lzw",
#      bg="white",
#      res=600)   
#PDF
pdf(file = "ROCs_all_5yr_291train.pdf", width = 8, height = 8)

par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt_GSig_PFS_PANCANCER_train$days_to_last_follow_up, status=rt_GSig_PFS_PANCANCER_train$vital_status, marker = rt_GSig_PFS_PANCANCER_train$riskScore_Gsig_PFS_PANCANCER_train, predict.time =5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1.05),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("Gene signature only"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1)

#绘制ROC曲线
roc=survivalROC(Stime=rt_Gleason_PFS_PANCANCER_train$days_to_last_follow_up, status=rt_Gleason_PFS_PANCANCER_train$vital_status, marker = rt_Gleason_PFS_PANCANCER_train$riskScore_GS_train, predict.time =5, method="KM")
aucText=c(aucText,paste0("Gleason score only"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1.05),col=rocCol[2],lwd = 2)

#绘制ROC曲线
roc=survivalROC(Stime=rt_GSigPlusGleason_PFS_PANCANCER_train$days_to_last_follow_up, status=rt_GSigPlusGleason_PFS_PANCANCER_train$vital_status, marker = rt_GSigPlusGleason_PFS_PANCANCER_train$riskScore_GSigPlusGS_train, predict.time =5, method="KM")
aucText=c(aucText,paste0("Gleason Score and gene signature"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1.05),col=rocCol[3],lwd = 2)

legend("bottomright", aucText,lwd=1,bty="n",col=rocCol)
dev.off()


## For testing cohort----

library(survivalROC)

## run together！！！
rocCol=c("#CCCCCC","#989898","#323232")
aucText=c()

#ROC curve
##TIFF
# tiff(file="ROCs_all_3yr.tiff",
#      width=18,
#      height=18,
#      units ="cm",
#      compression="lzw",
#      bg="white",
#      res=600)   
#PDF
pdf(file = "ROCs_all_5yr_125test.pdf", width = 8, height = 8)

par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt_GSig_PFS_PANCANCER_test$days_to_last_follow_up, status=rt_GSig_PFS_PANCANCER_test$vital_status, marker = rt_GSig_PFS_PANCANCER_test$riskScore_Gsig_PFS_PANCANCER_test, predict.time =5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1.05),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("Gene signature only"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1)

#绘制ROC曲线
roc=survivalROC(Stime=rt_Gleason_PFS_PANCANCER_test$days_to_last_follow_up, status=rt_Gleason_PFS_PANCANCER_test$vital_status, marker = rt_Gleason_PFS_PANCANCER_test$riskScore_GS_test, predict.time =5, method="KM")
aucText=c(aucText,paste0("Gleason score only"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1.05),col=rocCol[2],lwd = 2)

#绘制ROC曲线
roc=survivalROC(Stime=rt_GSigPlusGleason_PFS_PANCANCER_test$days_to_last_follow_up, status=rt_GSigPlusGleason_PFS_PANCANCER_test$vital_status, marker = rt_GSigPlusGleason_PFS_PANCANCER_test$riskScore_GSigPlusGS_test, predict.time =5, method="KM")
aucText=c(aucText,paste0("Gleason Score and gene signature"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1.05),col=rocCol[3],lwd = 2)

legend("bottomright", aucText,lwd=1,bty="n",col=rocCol)
dev.off()

##########################################################################################################################################################################



###################################################################################
#### 1/3/5-year ROC curves for Gene Signature only with 'survivalROC' (B&W) #######
##################################################################################
#ROC for gene signature only with survivalROC package; looks better and AUC slightly higher for the training set!!!

# For training----
library(survivalROC)

## run together！！！
rocCol=c("#CCCCCC","#989898","#323232")
aucText=c()

#ROC curve
##TIFF
# tiff(file="ROCs_all_3yr.tiff",
#      width=18,
#      height=18,
#      units ="cm",
#      compression="lzw",
#      bg="white",
#      res=600)   
#PDF
pdf(file = "ROCs_Gsig_allyrs_survivalROC_345train.pdf", width = 8, height = 8)

par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt_GSig_PFS_PANCANCER_train$days_to_last_follow_up, status=rt_GSig_PFS_PANCANCER_train$vital_status, marker = rt_GSig_PFS_PANCANCER_train$riskScore_Gsig_PFS_PANCANCER_train, predict.time =1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1.05),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("1-year PFS"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1)

#绘制ROC曲线
roc=survivalROC(Stime=rt_GSig_PFS_PANCANCER_train$days_to_last_follow_up, status=rt_GSig_PFS_PANCANCER_train$vital_status, marker = rt_GSig_PFS_PANCANCER_train$riskScore_Gsig_PFS_PANCANCER_train, predict.time =3, method="KM")
aucText=c(aucText,paste0("3-year PFS"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1.05),col=rocCol[2],lwd = 2)

#绘制ROC曲线
roc=survivalROC(Stime=rt_GSig_PFS_PANCANCER_train$days_to_last_follow_up, status=rt_GSig_PFS_PANCANCER_train$vital_status, marker = rt_GSig_PFS_PANCANCER_train$riskScore_Gsig_PFS_PANCANCER_train, predict.time =5, method="KM")
aucText=c(aucText,paste0("5-year PFS"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1.05),col=rocCol[3],lwd = 2)

legend("bottomright", aucText,lwd=1,bty="n",col=rocCol)
dev.off()


# For testing----
library(survivalROC)

## run together！！！
rocCol=c("#CCCCCC","#989898","#323232")
aucText=c()

#ROC curve
##TIFF
# tiff(file="ROCs_all_3yr.tiff",
#      width=18,
#      height=18,
#      units ="cm",
#      compression="lzw",
#      bg="white",
#      res=600)   
#PDF
pdf(file = "ROCs_Gsig_allyrs_survivalROC_147test.pdf", width = 8, height = 8)

par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt_GSig_PFS_PANCANCER_test$days_to_last_follow_up, status=rt_GSig_PFS_PANCANCER_test$vital_status, marker = rt_GSig_PFS_PANCANCER_test$riskScore_Gsig_PFS_PANCANCER_test, predict.time =1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1.05),col=rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText=c(aucText,paste0("1-year PFS"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1)

#绘制ROC曲线
roc=survivalROC(Stime=rt_GSig_PFS_PANCANCER_test$days_to_last_follow_up, status=rt_GSig_PFS_PANCANCER_test$vital_status, marker = rt_GSig_PFS_PANCANCER_test$riskScore_Gsig_PFS_PANCANCER_test, predict.time =3, method="KM")
aucText=c(aucText,paste0("3-year PFS"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1.05),col=rocCol[2],lwd = 2)

#绘制ROC曲线
roc=survivalROC(Stime=rt_GSig_PFS_PANCANCER_test$days_to_last_follow_up, status=rt_GSig_PFS_PANCANCER_test$vital_status, marker = rt_GSig_PFS_PANCANCER_test$riskScore_Gsig_PFS_PANCANCER_test, predict.time =5, method="KM")
aucText=c(aucText,paste0("5-year PFS"," (AUC=",sprintf("%.3f",roc$AUC),")"))
lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1.05),col=rocCol[3],lwd = 2)

legend("bottomright", aucText,lwd=1,bty="n",col=rocCol)
dev.off()





##################################################################
########################## Risk heatmap ##########################
##################################################################


## For training cohort----
rt_GSig_PFS_PANCANCER_train3=rt_GSig_PFS_PANCANCER_train[order(rt_GSig_PFS_PANCANCER_train$riskScore_Gsig_PFS_PANCANCER_train),]
rt1_train=rt_GSig_PFS_PANCANCER_train3[c(3:(ncol(rt_GSig_PFS_PANCANCER_train3)-2))]
rt1_train=t(rt1_train)

rt1_train=log2(rt1_train+1)
library(pheatmap)
annotation=data.frame(type=rt_GSig_PFS_PANCANCER_train3[,ncol(rt_GSig_PFS_PANCANCER_train3)])
rownames(annotation)=rownames(rt_GSig_PFS_PANCANCER_train3)

pdf(file="heatmap_train.pdf",width = 12,height = 5)
pheatmap(rt1_train, 
         annotation=annotation, 
         cluster_cols = FALSE,
         fontsize_row=11,
         fontsize_col=3,
         color = colorRampPalette(c("green", "black", "red"))(50) )
dev.off()

## For testing cohort----

rt_GSig_PFS_PANCANCER_test3=rt_GSig_PFS_PANCANCER_test[order(rt_GSig_PFS_PANCANCER_test$riskScore_Gsig_PFS_PANCANCER_test),]
rt1_test=rt_GSig_PFS_PANCANCER_test3[c(3:(ncol(rt_GSig_PFS_PANCANCER_test3)-2))]
rt1_test=t(rt1_test)

rt1_test=log2(rt1_test+1)
library(pheatmap)
annotation=data.frame(type=rt_GSig_PFS_PANCANCER_test3[,ncol(rt_GSig_PFS_PANCANCER_test3)])
rownames(annotation)=rownames(rt_GSig_PFS_PANCANCER_test3)

pdf(file="heatmap_test.pdf",width = 12,height = 5)
pheatmap(rt1_test, 
         annotation=annotation, 
         cluster_cols = FALSE,
         fontsize_row=11,
         fontsize_col=3,
         color = colorRampPalette(c("green", "black", "red"))(50) )
dev.off()



##################################################################
########################## Risk score ############################
##################################################################

## For training cohort----


#rt=rt[order(rt$riskScore),]

surv.cut_PRAD <- surv_cutpoint(
  rt_GSig_PFS_PANCANCER_train3,
  time = "days_to_last_follow_up",
  event = "vital_status",
  variables = "riskScore_Gsig_PFS_PANCANCER_train"
)



riskClass_train=rt_GSig_PFS_PANCANCER_train3[,"risk"]
lowLength=length(riskClass_train[riskClass_train=="low"])
highLength=length(riskClass_train[riskClass_train=="high"])
line=rt_GSig_PFS_PANCANCER_train3[,"riskScore_Gsig_PFS_PANCANCER_train"]
#line[line>10]=10
pdf(file="riskScore_train.pdf",width = 12,height = 5)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("green",lowLength),
           rep("red",highLength)))
legend("topleft",
       c("High risk", "Low risk"),
       lwd = 2,
       col = c("red","green"))
abline(h=summary(surv.cut_PRAD)[1],v=lowLength,lty=2)#median(rt$riskScore)
dev.off()




## For testing cohort----

#rt=rt[order(rt$riskScore),]
surv.cut_PRAD <- surv_cutpoint(
  rt_GSig_PFS_PANCANCER_test3,
  time = "days_to_last_follow_up",
  event = "vital_status",
  variables = "riskScore_Gsig_PFS_PANCANCER_test"
)



riskClass_test=rt_GSig_PFS_PANCANCER_test3[,"risk"]
lowLength=length(riskClass_test[riskClass_test=="low"])
highLength=length(riskClass_test[riskClass_test=="high"])
line=rt_GSig_PFS_PANCANCER_test3[,"riskScore_Gsig_PFS_PANCANCER_test"]
#line[line>10]=10
pdf(file="riskScore_test.pdf",width = 12,height = 5)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("green",lowLength),
           rep("red",highLength)))
legend("topleft",
       c("High risk", "Low risk"),
       lwd = 2,
       col = c("red","green"))
abline(h=summary(surv.cut_PRAD)[1],v=lowLength,lty=2)#median(rt$riskScore)
dev.off()



##################################################################
########################## Risk status ###########################
##################################################################

## For training cohort----
#rt=rt[order(rt$riskScore),]
riskClass_train=rt_GSig_PFS_PANCANCER_train3[,"risk"]
lowLength=length(riskClass_train[riskClass_train=="low"])
highLength=length(riskClass_train[riskClass_train=="high"])
color=as.vector(rt_GSig_PFS_PANCANCER_train3$vital_status)
color[color==1]="red"
color[color==0]="green"
pdf(file="survStat_train.pdf",width = 12,height = 5)
plot(rt_GSig_PFS_PANCANCER_train3$days_to_last_follow_up,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)

legend("topleft",
       c("Progression", "Not progression"),
       #lwd = 2,
       pch=16,
       col = c("red","green"))

abline(v=lowLength,lty=2)
dev.off()



## For testing cohort----


#rt=rt[order(rt$riskScore),]
riskClass_test=rt_GSig_PFS_PANCANCER_test3[,"risk"]
lowLength=length(riskClass_test[riskClass_test=="low"])
highLength=length(riskClass_test[riskClass_test=="high"])
color=as.vector(rt_GSig_PFS_PANCANCER_test3$vital_status)
color[color==1]="red"
color[color==0]="green"
pdf(file="survStat_test.pdf",width = 12,height = 5)
plot(rt_GSig_PFS_PANCANCER_test3$days_to_last_follow_up,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
legend("topleft",
       c("Progression", "Not progression"),
       #lwd = 2,
       pch=16,
       col = c("red","green"))

abline(v=lowLength,lty=2)
dev.off()

