##################################################################
############################# data import ########################
##################################################################
# load finalised ExpressionSet of 22165 genes
lnames_ExpressionSet <- load('/Users/zhuofanmou/All documents/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Gene identification/data_needed/palmieri_final_PC.RData')
lnames_ExpressionSet

##################################################################
################### Pairwise DGE analysis with 'Limma' #############
##################################################################

Biobase::pData(palmieri_final_PC)$Sample.ID_for_DGE = c('10','10','12','12','13','13','14','14','15_2','15_2','17','17','7','7','8','8','9','9')

individual_PC = as.character(Biobase::pData(palmieri_final_PC)$Sample.ID_for_DGE) #label the tissues with their corresponding patient numbers
tumor = Biobase::pData(palmieri_final_PC)$Cancer.status  #group_list
tumor = factor(tumor)
individual_PC = factor(individual_PC)
design_patasvar_PC <- model.matrix(~0+tumor+individual_PC) # design matrix
design_patasvar_PC <- design_patasvar_PC[ , c (2,1,3:10)] # re-arrange the columns
colnames(design_patasvar_PC)[1:2] <- c("T", "B")
rownames(design_patasvar_PC) <- individual_PC 
head(design_patasvar_PC)

contrast_matrix_PC <- makeContrasts( T-B, levels = design_patasvar_PC) #contrast matrix
contrast_matrix_PC


palmieri_fit_PC <- eBayes(contrasts.fit(lmFit(palmieri_final_PC,
                                              design = design_patasvar_PC),
                                        contrast_matrix_PC))
# results <- decideTests (palmieri_fit_PC)
# summary (results)

table_PC <- topTable(palmieri_fit_PC, number = Inf, adjust.method = "BH")  #results
#head(table_PC)
table_PC[1:6,5:10]

table(decideTests(palmieri_fit_PC, adjust.method = "BH", p.value = 0.1, lfc = log2(1.5))) # table of gene classes; '-1' = downregulated, '1' = upregulated, '0' = not

table_PC = na.omit(table_PC)    #！ This is the result we want for all 22165 genes!!!
openxlsx::write.xlsx(table_PC, file = "Limma_results_all.xlsx", append = FALSE) # use 'openxlsx' light memory and no need for java
#write.csv(table_PC,"limma_table_PC_results.csv",quote = F) # export results after DGE analysis for ALL genes
save(table_PC, file = "table_PC.RData") # saving limma results for all 22165 genes













