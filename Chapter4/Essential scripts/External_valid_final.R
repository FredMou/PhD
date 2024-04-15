

setwd("/Users/zhuofanmou/Downloads/survival_analysis")
getwd()


Signature_events_PSI_GSE54460 <- read.delim("~/Documents/PhD Exeter/ClariomD/alternative_splicing/R_workspace_v2/survival_analysis/Signature_events_PSI_GSE54460.txt")
View(Signature_events_PSI_GSE54460)
Signature_events_PSI_GSE54460 <- t(Signature_events_PSI_GSE54460)
colnames(Signature_events_PSI_GSE54460) <- Signature_events_PSI_GSE54460[1 , ]
Signature_events_PSI_GSE54460 <- Signature_events_PSI_GSE54460[ -1 , ]

table(is.na(Signature_events_PSI_GSE54460))

#Replace na values with 0 using is.na()
Signature_events_PSI_GSE54460[is.na(Signature_events_PSI_GSE54460)] = 0

#Display the dataframe
print(Signature_events_PSI_GSE54460)















################################################################################################################
########################################### External validation ################################################
################################################################################################################
GSE107299 <- load("isoform_PSI_NEBC_mle_QN_RMA.rda")
GSE107299
View(PSI)
PSI <- as.data.frame(PSI)
#PSI_final <- PSI[(rownames(PSI) == "TC06002344.hg_1") , ]
#PSI_final <- PSI[(rownames(PSI) == "TC19000275.hg_6" | rownames(PSI) == "TC19001901.hg_1" | rownames(PSI) == "TC07003299.hg_4" | rownames(PSI) == "TC03001359.hg_2" | rownames(PSI) == "TC19000458.hg_5") , ]

#PSI_final <- PSI[(rownames(PSI) == "TC19000275.hg_6" | rownames(PSI) == "TC07003299.hg_4" | rownames(PSI) == "TC19000458.hg_5") , ] # NOTE: the order of the colnames/probenames after extracting will not be in order, so later when assigning the colnames, we need to be careful!!!!!!! Need to double check 'PSI_final'


PSI_final <- PSI[(rownames(PSI) == "TC19000275.hg_6" | rownames(PSI) == "TC14000168.hg_10" | rownames(PSI) == "TC09001050.hg_4"  | rownames(PSI) == "TC07003299.hg_4" | rownames(PSI) == "TC03001359.hg_2" | rownames(PSI) == "TC19000458.hg_5") , ] # NOTE: the order of the colnames/probenames after extracting will not be in order, so later when assigning the colnames, we need to be careful!!!!!!! Need to double check 'PSI_final'



#PSI_final <- PSI[(rownames(PSI) == "TC19000275.hg_6" | rownames(PSI) == "TC19001901.hg_1" | rownames(PSI) == "TC07003299.hg_4" | rownames(PSI) == "TC03001359.hg_2" | rownames(PSI) == "TC19000458.hg_5") , ]



#PSI_final <- PSI[(rownames(PSI) == "TC19000527.hg_2" | rownames(PSI) == "TC06002344.hg_1" | rownames(PSI) == "TC17000565.hg_9") , ]

PSI_final <- t(PSI_final) 
View(PSI_final)





GSE107299_cli <- readRDS("CPC-Gene_eSet.RDS")
GSE107299_cli
GSE107299_cli <- pData(GSE107299_cli)
GSE107299_cli <- GSE107299_cli[ GSE107299_cli$batch == "Batch2", ]
GSE107299_cli <- GSE107299_cli[ , c(26,27)]
GSE107299_cli$time_to_bcr <- GSE107299_cli$time_to_bcr/12 
View(GSE107299_cli)

## Combine the data 
GSE107299_final <- cbind(GSE107299_cli,PSI_final) 
#GSE107299_final <- GSE107299_final[ , -3]
#colnames(GSE107299_final)[1:7] <- c("futime" , "fustat" , "ALS2CL.64461.RI","CYP3A5.80711.RI","CYP4F12.48110.RI","FXYD3.49039.RI","ZNF154.52329.RI")

#colnames(GSE107299_final)[1:6] <- c("futime" , "fustat" , "CYP3A5.80711.RI", "NFATC4.26991.RI", "CYP4F12.48110.RI","FXYD3.49039.RI")

colnames(GSE107299_final)[1:8] <- c("futime" , "fustat" , "ALS2CL.64461.RI", "CYP3A5.80711.RI", "PIGO.86233.RI", "NFATC4.26991.RI", "CYP4F12.48110.RI","FXYD3.49039.RI")

#colnames(GSE107299_final)[1:5] <- c("futime" , "fustat" , "CYP3A5.80711.RI", "CYP4F12.48110.RI","FXYD3.49039.RI")

GSE107299_final <- na.omit(GSE107299_final)
View(GSE107299_final)


#GSE107299_final <- GSE107299_final[ , c(1,2,5,4,3,6)]

######################################### LASSO-COX direct model prediction
# Ref: 
# https://github.com/rli012/PCaSignatures/blob/master/Evaluation_Survival_Analysis_Intra_Dataset.R
# https://github.com/rli012/PCaSignatures/blob/master/Evaluation_Survival_Analysis_Inter_Dataset.R
#riskScore <- predict(cvfit, s=cvfit$lambda.min, newx = x_test, type="link") # response: reltive risk; link: linear prediction



# multiCox_GSE107299 <- coxph(Surv(futime, fustat) ~ ., data = GSE107299_final)#coxph(Surv(futime, fustat) ~ `PDCD2|78501|RI`, data = rt)#coxph(Surv(futime, fustat) ~ ., data = rt)
# multiCox_GSE107299



# # calculate the risk score manually from the external dataset
#riskScore <-  1.293*GSE107299_final$CYP4F12.48110.RI + 2.09*GSE107299_final$ZNF154.52329.RI + 0.9199*GSE107299_final$CYP3A5.80711.RI + 1.206*GSE107299_final$ALS2CL.64461.RI + 11.69*GSE107299_final$FXYD3.49039.RI

#riskScore <-  1.882*GSE107299_final$CYP4F12.48110.RI + 1.919*GSE107299_final$NFATC4.26991.RI + 1.187*GSE107299_final$CYP3A5.80711.RI + 12.65*GSE107299_final$FXYD3.49039.RI


riskScore <- 1.374*GSE107299_final$CYP4F12.48110.RI + 1.844*GSE107299_final$NFATC4.26991.RI + (-4.517)*GSE107299_final$PIGO.86233.RI + 0.9702*GSE107299_final$CYP3A5.80711.RI + 1.013*GSE107299_final$ALS2CL.64461.RI + 15.59*GSE107299_final$FXYD3.49039.RI


#riskScore<- predict(multiCox, type="lp", newdata=GSE107299_final, se.fit = FALSE) # or use the model developed directly, the ROC analysis will be the the same as the riskScore calculated above

#riskScore <- 3.719*GSE107299_final$PDCD2.78501.RI
rt <- GSE107299_final
#riskScore<- predict(multiCox, type="lp", newdata=rt, se.fit = FALSE)
#View(as.data.frame(riskScore))


outCol=colnames(GSE107299_final)
rt <- GSE107299_final
rt$risk <- riskScore
View(rt)
class(rt)

library(maxstat)

opt_cut <-  maxstat.test(Surv(futime,fustat) ~ risk, data = rt, smethod = "LogRank",pmethod = "Lau92")#, pmethod = "Lau92")
opt_cut
cutpoint<- opt_cut$estimate
cutpoint

# # Below use the best cut-off as riskscore cutoff:
# training.risk.score <- riskScore[ , 1]
# class(training.risk.score)
# training.risk.group <- training.risk.score > cutpoint
# table(training.risk.group)


risk=as.vector(ifelse(riskScore>cutpoint, "high", "low"))
table(risk)
outTab=cbind(rt[,outCol], riskScore=as.vector(riskScore), risk)
write.table(cbind(id=rownames(outTab),outTab), file="risk_GSE107299.txt", sep="\t", quote=F, row.names=F)

risk=read.table("risk_GSE107299.txt", header=T, sep="\t", check.names=F, row.names=1)
View(risk)

diff <- survdiff(Surv(futime,fustat) ~ risk, data = risk)
pValue <- 1-pchisq(diff$chisq, df=1)
pValue <- signif(pValue, 4)
pValue <- format(pValue, scientific=FALSE)



fit_PRAD <- survfit(Surv(futime,fustat) ~ risk,
                    data = risk)

pdf(file = "survival_optimalBased_signatureValid_GSE107299.pdf", width = 8, height = 8)

surPlot_bestcut=ggsurvplot(fit_PRAD,
                           data=risk, #这个可要可不要
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






#引用包
library(survival)
library(survminer)
library(timeROC)
riskFile="risk_GSE107299.txt"         #风险输入文件
#cliFile="clinical_categorised.txt"      #临床数据文件
#setwd("C:\\biowolf\\immAS\\16.ROC")     #修改工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
View(risk)
risk=risk[,c("futime", "fustat", "riskScore")]

#读取临床数据文件
# cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
# View(cli)
# str(cli)
# #合并数据
# samSample=intersect(row.names(risk), row.names(cli))
# risk1=risk[samSample,,drop=F]
# cli=cli[samSample,,drop=F]
rt=risk
View(rt)
median(rt$futime)
# NEW BIT ADDING FOR VARIABLE TRANSFORMATION (BUT DOES NOT MATTER FOR 2-LEVEL CATEGORICAL VARIABLE I THINK!!! ANSWER: YES!!!TRUE FOR BOTH UNI&MULTI)
# It needs to be number 1 or 0 or more or continuous (for sure) for ROC curve analysis
# rt[rt$age == '<60',]$age <- 0
# rt[rt$age == '>=60',]$age <- 1

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

#定义颜色
bioCol=rainbow(ncol(rt)-1, s=0.9, v=0.9)


######绘制 N 年的ROC曲线######
ROC_rt=timeROC(T=risk$futime,delta=risk$fustat,
               marker=risk$riskScore,cause=1,
               weighting='cox',
               times=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),ROC=TRUE)
ROC_rt$AUC
####### For 3 and 5-year ROC only
pdf(file="ROC_358yrs_GSE107299.pdf", width=5.5, height=5.5)
plot(ROC_rt,time=3,col="#27347F",title=FALSE,lwd=2)
plot(ROC_rt,time=5,col="#AB1815",add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=8,col="#008000",add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[3])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[5])),
         paste0('AUC at 8 years: ',sprintf("%.03f",ROC_rt$AUC[8]))),
       
       col=c("#27347F","#AB1815","#008000"), lwd=2, bty = 'n')
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
