

# setwd("/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/KM curves for clinical groups")
# getwd()

##################################################################
####################### MulCox construction ######################
##################################################################
lnames_trainset <- load('/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Model construction/R_data_needed/train.rda')
lnames_trainset
lnames_testset <- load('/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Model construction/R_data_needed/test.rda')
lnames_testset
# Inputting the selected gene signature from LASSO----
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
1-with(train, rcorr.cens(riskScore_Gsig_PFS_PANCANCER_train, Surv(days_to_last_follow_up,vital_status)))[1] 
1-with(test, rcorr.cens(riskScore_Gsig_PFS_PANCANCER_test, Surv(days_to_last_follow_up,vital_status)))[1] 

# optiml cut-off to classify riskscore on the Train and Test set now!----
## I used NA remained clinical data!!! And need to run mannually/fsubfactor-by-subfctor
# For training set
lnames_train <- load('/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Model construction/R_data_needed/df_clinic.genemarkers.COX_train.rda')
lnames_train

KM_train <-df_clinic.genemarkers.COX_train
KM_train$risk <- riskScore_Gsig_PFS_PANCANCER_train
View(KM_train)

KM_train_gleason <- subset(KM_train,KM_train$gleason_score == "≥8") # CHANGE:<8
View(KM_train_gleason)

KM_train_T <- subset(KM_train,KM_train$ajcc_pathologic_t == "T3-4") #CHANGE:T2
View(KM_train_T)

KM_train_N <- subset(KM_train,KM_train$ajcc_pathologic_n == "N1") #CHANGE:N0
View(KM_train_N)

# For testing set
lnames_test <- load('/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Model construction/R_data_needed/df_clinic.genemarkers.COX_test.rda')
lnames_test

KM_test <-df_clinic.genemarkers.COX_test
KM_test$risk <- riskScore_Gsig_PFS_PANCANCER_test
View(KM_test)

KM_test_gleason <- subset(KM_test,KM_test$gleason_score == "≥8") #CHANGE:<8
View(KM_test_gleason)

KM_test_T <- subset(KM_test,KM_test$ajcc_pathologic_t == "T3-4")  #CHANGE:T2
View(KM_test_T)

KM_test_N <- subset(KM_test,KM_test$ajcc_pathologic_n == "N1")  #CHANGE:N0
View(KM_test_N)



# Optimal cut-ff approach for survival analysis of calculated riskScore ----
surv.cut_PRAD <- surv_cutpoint(
  KM_train_gleason, # CHANGE
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

#KM analysis and survival curve based on the high/low risk
fit_PRAD <- survfit(Surv(days_to_last_follow_up,vital_status) ~ risk,
                    data = surv.cat_PRAD)

pdf(file = "survival_MulCox_KMplot_GSig_Gleason_GE8_train.pdf", width = 8, height = 8) #CHANGE

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

