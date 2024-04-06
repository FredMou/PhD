lnames <- load('/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/data_prad.rda')
lnames






library(edgeR)
exp_edgeR <- read.table('/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/datExp.counts.prad.txt', header=T, sep="\t", check.names=F)
dim(exp_edgeR)
exp_edgeR=as.data.frame(exp_edgeR)
dim(exp_edgeR)
anno = read.csv("/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/Annotation.csv",header = T)
exp_edgeR$Gene = anno$Gene.Symbol[match(rownames(exp_edgeR),anno$ID)] # gene symbol added to the last (552) column
dim(exp_edgeR)
table(is.na(exp_edgeR$Gene)) # 4250 missing/NA genes! 4205(NA) vs 52397(annotated)
View(exp_edgeR$Gene)
exp_edgeR <- na.omit(exp_edgeR)
dim(exp_edgeR)
exp_edgeR=as.matrix(exp_edgeR)

rownames(exp_edgeR)=exp_edgeR[,552]
View(head(exp_edgeR))
dim(exp_edgeR)
exp_edgeR=exp_edgeR[,-552]
dimnames=list(rownames(exp_edgeR),colnames(exp_edgeR))
dim(exp_edgeR) #52397*551
View(exp_edgeR)
class(exp_edgeR)

table(duplicated(rownames(exp_edgeR))) #1388个重复的基因名
any(duplicated(rownames(exp_edgeR)))


data_edgeR=matrix(as.numeric(as.matrix(exp_edgeR)), nrow=nrow(exp_edgeR), dimnames=dimnames)
class(data_edgeR)


data_edgeR=avereps(data_edgeR) # taking average if a gene is duplicated
dim(data_edgeR) #51009*551
View(data_edgeR)
#data_edgeR <- data_edgeR[ rowMeans(data_edgeR)>1, ] #29821*551


any(duplicated(rownames(data_edgeR))) # No duplicated genes now!

exp_edgeR <- t(data_edgeR) #这里转置一下是为了后面的group， 没别的用了因为后面还是会转置的
dim(exp_edgeR)  #551*51009

#Clinic info all 551 samples----
clinic_info_prad <- sample.clinic_prad
dim(clinic_info_prad)
group <- clinic_info_prad[ , c('barcode','shortLetterCode')]
group[group$shortLetterCode == 'TP',]$shortLetterCode <- 'tumor' 
group[group$shortLetterCode == 'TM',]$shortLetterCode <- 'tumor' 
group[group$shortLetterCode == 'NT',]$shortLetterCode <- 'normal' 

# check for sample order match
table(rownames(exp_edgeR)==rownames(group)) #一开始就是顺序对上了！

group <- group$shortLetterCode
class(group)
View(group)

#edgeR normalisation
#https://www.bilibili.com/video/BV1jo4y1D74s?from=search&seid=17345839985870354983
#https://www.bilibili.com/video/BV18E411Q73b?p=8
datNorm_edgeR <- DGEList(counts = data_edgeR, group = group)
keep <- rowSums(cpm(datNorm_edgeR)>1) >= 2
datNorm_edgeR <- datNorm_edgeR[keep, , keep.lib.sizes = FALSE]
dim(datNorm_edgeR)
datNorm_edgeR <- calcNormFactors(datNorm_edgeR)
datNorm_edgeR <- estimateCommonDisp(datNorm_edgeR)
datNorm_edgeR <- estimateTagwiseDisp(datNorm_edgeR)
datNorm_all_edgeR <- datNorm_edgeR$pseudo.counts
#datNorm_all_edgeR <- rbind(id = colnames(datNorm_all_edgeR), datNorm_all_edgeR) #没必要，可能对输出友好，就多加了一行样本名且以id为名
View(head(datNorm_all_edgeR))
dim(datNorm_all_edgeR)


datNorm_all_edgeR <- t(datNorm_all_edgeR)
save(datNorm_all_edgeR, file = "ExpDat_norm_edgeR.rda") #Store normalised data


datNorm_all_log2_edgeR <- log2(datNorm_all_edgeR[ , 1:ncol(datNorm_all_edgeR)]+1)
View(head(datNorm_all_log2_edgeR))
dim(datNorm_all_log2_edgeR)
save(datNorm_all_log2_edgeR, file = "ExpDat_norm_log2_edgeR.rda") #Store log2-transformed normalised data

# five_gene_edgeR <-datNorm_all_log2_edgeR[ , c('HNRNPL', 'ELAVL1', 'PCBP1', 'PCBP2', 'PABPN1','PTPRF','MRPL24','DANCR', 'MYC', 'TRPM4')]#exp_fkpm[,c('MRPL24','MYC','HNRNPL','ELAVL1','IMPDH2')]
# dim(five_gene_edgeR)
# 
# 
# #Clinic info tumour samples----
# #final complete clinical info with survival time added for some dead samples
# clinic_info_Cox <- read.table("/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/clinical.finalCox.txt.txt", header = T, sep = "\t") 
# 
# dim(clinic_info_Cox)
# View(clinic_info_Cox)
# 
# 
# 
# pheno_edgeR <- clinic_info_Cox
# 
# dim(pheno_edgeR)
# View(pheno_edgeR)
# pheno_edgeR <- pheno_edgeR [ -(pheno_edgeR$shortLetterCode=='TM'),]
# pheno_edgeR <- pheno_edgeR[c('barcode','days_to_last_follow_up','vital_status')]#对于edgeR， 我们要留‘barcode’名
# rownames(pheno_edgeR) <- pheno_edgeR$barcode
# dim(pheno_edgeR)
# 
# #合并----
# #check rowname match
# sameSample_edgeR=intersect(row.names(pheno_edgeR),row.names(five_gene_edgeR)) #以'pheno_survival_genesigOnly' 的行名来找因为它在前
# length(sameSample_edgeR)
# head(sameSample_edgeR)
# 
# pheno_edgeR2=pheno_edgeR[sameSample_edgeR,]
# dim(pheno_edgeR2)
# head(rownames(pheno_edgeR2))
# 
# five_gene_edgeR2=five_gene_edgeR[sameSample_edgeR,]
# dim(five_gene_edgeR2)
# head(rownames(five_gene_edgeR2))
# 
# #Check if row names matched
# table(row.names(five_gene_edgeR2) == row.names(pheno_edgeR2))
# 
# 
# 
# 
# # data input prep of selected genes for survival analysis----
# CombindDat_model_edgeR=cbind(pheno_edgeR2,five_gene_edgeR2)
# dim(CombindDat_model_edgeR)
# View(CombindDat_model_edgeR)
# CombindDat_model_edgeR <- CombindDat_model_edgeR[ , -1]
# CombindDat_model_edgeR[CombindDat_model_edgeR$vital_status == 'Alive',]$vital_status <- 0 #Alive
# CombindDat_model_edgeR[CombindDat_model_edgeR$vital_status == 'Dead',]$vital_status <- 1 #Dead
# CombindDat_model_edgeR$vital_status[!CombindDat_model_edgeR$vital_status %in% c("0","1")] = NA
# class(CombindDat_model_edgeR$vital_status)#class(CombindDat_model$gleason_score)
# CombindDat_model_edgeR$vital_status <- as.numeric(CombindDat_model_edgeR$vital_status)#CombindDat_model$gleason_score <- as.factor(CombindDat_model$gleason_score)
# table(CombindDat_model_edgeR$vital_status)
# View(CombindDat_model_edgeR$vital_status)
# # colnames(pheno) = c("barcode","time","status")
# ## FINAL data frame of selected gene for survival:
# CombindDat_model_edgeR$days_to_last_follow_up = CombindDat_model_edgeR$days_to_last_follow_up / 365 # /30 for month;
# class(CombindDat_model_edgeR$vital_status)
# class(CombindDat_model_edgeR$days_to_last_follow_up)
# dim(CombindDat_model_edgeR)
# View(CombindDat_model_edgeR)
# CombindDat_model_edgeR <- na.omit(CombindDat_model_edgeR) # 这里如果是用前面基因的话是不会有确实值的，所以跑不跑无所谓
# 
# ## Adding other cilinicpathological info:
# 
# 
# 
# 
# 
# 
# 
# # Multi Cox----
# class(CombindDat_model_edgeR)
# dim(CombindDat_model_edgeR)
# cox_Mul <- coxph(Surv(days_to_last_follow_up,vital_status) ~ MYC+ELAVL1  , data = CombindDat_model_edgeR) #age_at_index+ajcc_pathologic_t+ajcc_pathologic_n+gleason_score+MRPL24+MYC+HNRNPL+ELAVL1+IMPDH2#coxph(Surv(days_to_last_follow_up,vital_status) ~ ., data = CombindDat_model)
# #cox_Mul <- step(cox_Mul, direction = "both") # 几乎什么情况都是gleason 被选出; 上边一步中直接写出这步得到变量也可以的，结果一样
# riskScore_genesig_clinic <- predict(cox_Mul, type = "risk", newdata = CombindDat_model_edgeR)
# summary <- summary(cox_Mul)
# coxGene <- rownames(summary$coefficients)
# coxGene <- gsub("","",coxGene)
# outCol <- c("days_to_last_follow_up","vital_status",'RPLP0', 'TPI1', 'IMPDH2', 'MRPL24','APEX1','HNRNPL', 'ELAVL1', 'PCBP1', 'PCBP2', 'PABPN1','PTPRF', 'DANCR', 'MYC', 'TRPM4')
# View(outCol)
# #outCol <- "gleason_score"
# risk <- as.vector(ifelse(riskScore_genesig_clinic>median(riskScore_genesig_clinic), "high","low"))
# write.table(cbind(id = rownames(cbind(CombindDat_model_edgeR[ , outCol], riskScore_genesig_clinic, risk)), cbind(CombindDat_model_edgeR[ , outCol],riskScore_genesig_clinic,risk)),
#             file = "risk.txt",
#             sep = "\t",
#             quote = F,
#             row.names = F)
# 
# rt <- read.table("/Users/zhuofanmou/All documents/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/risk.txt", header = T, sep = "\t", check.names = F, row.names = 1) 
# View(rt)
# 
# 
# #ROC
# tiff(file = "MulCox_ROC_10yr_5gene_fpkm.tiff",
#      width = 30,
#      height = 30,
#      units = "cm",
#      compression = "lzw",
#      bg = "white",
#      res = 600)
# 
# par(oma = c(0.5,1,0,1), font.lab = 1.5, font.axis = 1.5)
# 
# roc <- survivalROC(Stime = rt$days_to_last_follow_up, status = rt$vital_status, marker = rt$riskScore_genesig_clinic,
#                    predict.time = 1, method = "KM")
# 
# plot(roc$FP, roc$TP, type = "l", xlim = c(0,1), ylim = c(0,1), col = 'black',
#      xlab = "False positive rate", ylab = "True positive rate",
#      main = paste("ROC curve with 5genes  (10 year","AUC = ",round(roc$AUC,3),")"),
#      lwd = 2, cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.2, font = 1.2)
# 
# abline(0,1)
# 
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# #single gene analysis using Kaplan-M (but as a batch)----
# picDir <- 'KM plot 14genes edgeR'
# dir.create(picDir)
# 
# library(survival)
# outTab <- data.frame()
# 
# #K-M plot for all 14 genes at once
# for (gene in colnames(CombindDat_model_edgeR[ ,3:ncol(CombindDat_model_edgeR) ])){
#   a = CombindDat_model_edgeR[ , gene] < median(CombindDat_model_edgeR[ , gene])
#   diff = survdiff(Surv(days_to_last_follow_up,vital_status) ~ a, data = CombindDat_model_edgeR)
#   pValue = 1-pchisq(diff$chisq, df=1)
#   outTab = rbind(outTab, cbind(gene=gene, pvalue=pValue))
#   pValue=signif(pValue,4)
#   pValue=format(pValue,scientific=TRUE)
#   
#   fit <- survfit(Surv(days_to_last_follow_up,vital_status) ~ a, data = CombindDat_model_edgeR)
#   summary(fit)
#   
#   tiff(file = paste(picDir, "/", gene, ".survival.tiff", sep = ""),
#        width = 14,
#        height = 14,
#        units = "cm",
#        compression = "lzw",
#        bg = "white",
#        res = 600)
#   
#   plot(fit, 
#        lwd = 2,
#        col = c("red","blue"),
#        xlab = "Time (year)",
#        mark.time = T,
#        ylab = "Survival rate",
#        main = paste("Overall survival curve (p=", pValue ,")", sep = ""))
#   legend("topright",
#          c(paste(gene," high expression",sep = ""),
#            paste(gene," low expression",sep = "")),
#          lwd = 2,
#          col = c("red","blue"))
#   dev.off()
# }
# 
# #write.table(outTab, file = "survival.xls", sep = "\t", row.names = F, quote = F)
# 
# write.table(outTab, file = "survival_km_14genes_edgeR.txt", sep = "\t", row.names = F, quote = F, col.names = T)
# 
# 
# 
# 
# 
# x <- CombindDat_model_edgeR[ ,3:ncol(CombindDat_model_edgeR) ]
# dim(x)
# View(x)
# class(x)
# y <- CombindDat_model_edgeR[ , 2]
# class(y)
# length(y)
# dim(y)
# library(glmnet)
# 
# 
# #调优参数
# set.seed(1000) 
# #不设置的话，每次运行时后面的结果（选出的基因）会不断变化。这也就是说，如果文献中使用了Lasso回归，是没有办法复现结果的。
# cv_fit <- cv.glmnet(x=x, y=y) #cv.glmnet()就是一个调优参数的过程
# plot(cv_fit)
# 
# 
# #系数图
# fit <- glmnet(x=x, y=y) #glmnet是构建模型的
# plot(fit,xvar = "lambda")
# 
# 
# model_lasso_min <- glmnet(x=x, y=y,lambda=cv_fit$lambda.min) # 得到的模型是一个列表
# model_lasso_1se <- glmnet(x=x, y=y,lambda=cv_fit$lambda.1se) # 得到的模型是一个列表
# View(model_lasso_1se) #从中选一个查看一下格式
# 
# 
# head(model_lasso_min$beta,20)
# choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0] #as.numeric不等于0的基因就是有系数的也就是被选中了的基因。
# choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]
# length(choose_gene_min)
# View(choose_gene_min)
# # [1] 35 ##选到了35个基因
# length(choose_gene_1se)
# # [1] 11 ##选到了11个基因
# save(choose_gene_min,file = "lasso_choose_gene_min.Rdata")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #Univariate COX----
# outTab_UniCOX <- data.frame()
# library(survival)
# table(CombindDat_model_edgeR$vital_status)
# 
# for(i in colnames(CombindDat_model_edgeR[, 3:ncol(CombindDat_model_edgeR)])){
#   cox <- coxph(Surv(days_to_last_follow_up,vital_status) ~ CombindDat_model_edgeR[ , i], data = CombindDat_model_edgeR)
#   coxSummary = summary(cox)
#   outTab_UniCOX <- rbind(outTab_UniCOX, cbind(gene=i, HR=coxSummary$coefficients[,"exp(coef)"], 
#                                               z=coxSummary$coefficients[,"z"], 
#                                               pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
# }
# 
# write.table(outTab_UniCOX, file = "UnivariateCox_14genes.txt", sep = "\t", row.names = F, quote = F,col.names = T) # none significant; if some significant, screen the gene names manually
# 
# 
# 
# # Multi Cox
# class(CombindDat_model_edgeR)
# dim(CombindDat_model_edgeR)
# cox_Mul <- coxph(Surv(days_to_last_follow_up,vital_status) ~ . , data = CombindDat_model_edgeR) #age_at_index+ajcc_pathologic_t+ajcc_pathologic_n+gleason_score+MRPL24+MYC+HNRNPL+ELAVL1+IMPDH2#coxph(Surv(days_to_last_follow_up,vital_status) ~ ., data = CombindDat_model)
# cox_Mul <- step(cox_Mul, direction = "both") # 几乎什么情况都是gleason 被选出; 上边一步中直接写出这步得到变量也可以的，结果一样
# riskScore_genesig_clinic <- predict(cox_Mul, type = "risk", newdata = CombindDat_model_edgeR)
# summary <- summary(cox_Mul)
# coxGene <- rownames(summary$coefficients)
# coxGene <- gsub("","",coxGene)
# outCol <- c("days_to_last_follow_up","vital_status",'RPLP0', 'TPI1', 'IMPDH2', 'MRPL24','APEX1','HNRNPL', 'ELAVL1', 'PCBP1', 'PCBP2', 'PABPN1','PTPRF', 'DANCR', 'MYC', 'TRPM4')
# View(outCol)
# #outCol <- "gleason_score"
# risk <- as.vector(ifelse(riskScore_genesig_clinic>median(riskScore_genesig_clinic), "high","low"))
# write.table(cbind(id = rownames(cbind(CombindDat_model_edgeR[ , outCol], riskScore_genesig_clinic, risk)), cbind(CombindDat_model_edgeR[ , outCol],riskScore_genesig_clinic,risk)),
#             file = "risk.txt",
#             sep = "\t",
#             quote = F,
#             row.names = F)
# 
# rt <- read.table("/Users/zhuofanmou/All documents/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/risk.txt", header = T, sep = "\t", check.names = F, row.names = 1) 
# View(rt)
# 
# 
# #ROC
# library(survivalROC)
# tiff(file = "MulCox_ROC_5yr_14gene.tiff",
#      width = 30,
#      height = 30,
#      units = "cm",
#      compression = "lzw",
#      bg = "white",
#      res = 600)
# 
# par(oma = c(0.5,1,0,1), font.lab = 1.5, font.axis = 1.5)
# 
# roc <- survivalROC(Stime = rt$days_to_last_follow_up, status = rt$vital_status, marker = rt$riskScore_genesig_clinic,
#                    predict.time = 5, method = "KM")
# 
# plot(roc$FP, roc$TP, type = "l", xlim = c(0,1), ylim = c(0,1), col = 'black',
#      xlab = "False positive rate", ylab = "True positive rate",
#      main = paste("ROC curve with 14genes  (5 year","AUC = ",round(roc$AUC,3),")"),
#      lwd = 2, cex.main = 1.3, cex.lab = 1.2, cex.axis = 1.2, font = 1.2)
# 
# abline(0,1)
# 
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
