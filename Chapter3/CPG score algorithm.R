
setwd('/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v15/CPG score algorithm')
getwd()

# CPG algorithm construction from the curated power query analysis----

##For TCGA-PRAD training 
CPG_train = openxlsx::read.xlsx("TCGA.PRAD.training.CPG.final.xlsx")
rownames(CPG_train) <- CPG_train$Column1
CPG_train <- CPG_train[ , -1]
View(CPG_train)


CPG_train$score12 <- CPG_train$CPG1 + CPG_train$CPG2 # Numeric
table(CPG_train$score12)
CPG_train$score12[CPG_train$score12 == 3] <- 2  

CPG_train$score123 <- CPG_train$score12 + CPG_train$CPG3
table(CPG_train$score123)
CPG_train$score123[CPG_train$score123 == 4] <- 3
CPG_train$score123[CPG_train$score123 == 5] <- 3

CPG_train$score1234 <- CPG_train$score123 + CPG_train$CPG4
table(CPG_train$score1234)
CPG_train$score1234[CPG_train$score1234 == 5] <- 4
CPG_train$score1234[CPG_train$score1234 == 6] <- 4
CPG_train$score1234[CPG_train$score1234 == 7] <- 4


CPG_train$score12345 <- CPG_train$score1234 + CPG_train$CPG5
table(CPG_train$score12345)
CPG_train$score12345[CPG_train$score12345 == 6] <- 5
CPG_train$score12345[CPG_train$score12345 == 7] <- 5
CPG_train$score12345[CPG_train$score12345 == 8] <- 5
CPG_train$score12345[CPG_train$score12345 == 9] <- 5
table(CPG_train$score12345)  # final CPG score

# extracting variables for Cox model construction
df_CPG_COX_train <- CPG_train[ , c(1,2,8,12,13,14,37)]
View(df_CPG_COX_train)
colnames(df_CPG_COX_train)[7] <- "CPG_score"   # final training dataframe

# Create CPG risk table for supplementary 
supCPG_train <- df_CPG_COX_train
supCPG_train$SampleID <- rownames(supCPG_train) 
supCPG_train <- supCPG_train[ , c(8,7)]
View(supCPG_train)

openxlsx::write.xlsx(supCPG_train, file = "CPG_train_supplementary.xlsx", append = FALSE) # use 'openxlsx' light memory and no need for java
#write.csv(table_PC,"limma_table_PC_results.csv",quote = F) # export results after DGE analysis for ALL genes
#save(table_PC, file = "table_PC.RData") # saving limma results for all 22165 genes






##For TCGA-PRAD testing 
CPG_test = openxlsx::read.xlsx("TCGA.PRAD.testing.CPG.final.xlsx")
rownames(CPG_test) <- CPG_test$Column1
CPG_test <- CPG_test[ , -1]
View(CPG_test)


CPG_test$score12 <- CPG_test$CPG1 + CPG_test$CPG2 # Numeric
table(CPG_test$score12)
CPG_test$score12[CPG_test$score12 == 3] <- 2  

CPG_test$score123 <- CPG_test$score12 + CPG_test$CPG3
table(CPG_test$score123)
CPG_test$score123[CPG_test$score123 == 4] <- 3
CPG_test$score123[CPG_test$score123 == 5] <- 3

CPG_test$score1234 <- CPG_test$score123 + CPG_test$CPG4
table(CPG_test$score1234)
CPG_test$score1234[CPG_test$score1234 == 5] <- 4
CPG_test$score1234[CPG_test$score1234 == 6] <- 4
CPG_test$score1234[CPG_test$score1234 == 7] <- 4


CPG_test$score12345 <- CPG_test$score1234 + CPG_test$CPG5
table(CPG_test$score12345)
CPG_test$score12345[CPG_test$score12345 == 6] <- 5
CPG_test$score12345[CPG_test$score12345 == 7] <- 5
CPG_test$score12345[CPG_test$score12345 == 8] <- 5
CPG_test$score12345[CPG_test$score12345 == 9] <- 5
table(CPG_test$score12345) # final CPG score


# extracting variables for Cox model construction
df_CPG_COX_test <- CPG_test[ , c(1,2,8,12,13,14,37)]
View(df_CPG_COX_test)
colnames(df_CPG_COX_test)[7] <- "CPG_score" # final testing dataframe



# Create CPG risk table for supplementary 
supCPG_test <- df_CPG_COX_test
supCPG_test$SampleID <- rownames(supCPG_test) 
supCPG_test <- supCPG_test[ , c(8,7)]
View(supCPG_test)

openxlsx::write.xlsx(supCPG_test, file = "CPG_test_supplementary.xlsx", append = FALSE) # use 'openxlsx' light memory and no need for java



