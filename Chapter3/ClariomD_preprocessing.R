##################################################################
######################## working directory #######################
##################################################################
setwd('/Users/zhuofanmou/All documents/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R for 12.1/workspace')#Set this to your preferred working directory
##################################################################

##################################################################
######################## data import #############################
##################################################################
#CEL file----
CELs_PC <- list.celfiles("/Users/zhuofanmou/All documents/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Gene identification/data_needed/ClariomD_CEL", full.names = TRUE) # change this path to where the raw ClariomD CEL files are
raw_data_PC <- oligo::read.celfiles(CELs_PC) # read in the CEL files
stopifnot(validObject(raw_data_PC))




##################################################################
################## adding sample information #####################
##################################################################

Biobase::pData(raw_data_PC)$Sample.ID <- 
  c('10A1','10E5','12A1','12E5','13A1','13E5','14A1','14E5','15A1_2','15A1','15D4_2','15D4','17A1','17B2','7A1','7K11','8A1','8G7','904','9A1') # adding patients label
Biobase::pData(raw_data_PC)$Cancer.status <- 
  c('T','B','T','B','T','B','T','B','T','T','B','B','B','T','B','T','B','T','B','T') # adding Cancer status: T=maliganant and B=benign
head(exprs(raw_data_PC))
head(pData(raw_data_PC))
raw_data_15rem_PC = raw_data_PC[ , -c(10,12)] # removing duplicated patient 15 first pair
raw_data_PC = raw_data_15rem_PC
raw_data_PC


##################################################################
############ Robust multichip average (RMA) algorithm ############
##################################################################
# RMA algorithm on the raw data----
palmieri_eset_norm_PC <- oligo::rma(raw_data_PC, target = "core") # target='core' for transcript/gene level
palmieri_eset_norm_PC 

##################################################################
######################## Quality Check ###########################
##################################################################
# MA plot----
## on log2 raw data (manually transformed)
exprs_raw_PC = oligo::pm(raw_data_PC)
exprs_raw_log2_PC = log2(exprs_raw_PC)
affy::mva.pairs(exprs_raw_log2_PC[,c(4,8,12,18)],plot.method="normal",  main = "MAplot for Prostate log2 raw data") # Note that we only consider the perfect matched probes 

## on RMA normalised data
exprs_RMA_PC = exprs(palmieri_eset_norm_PC)
affy::mva.pairs(exprs_RMA_PC[,c(4,8,12,18)],plot.method="normal",  main = "MAplot for prostate data after RMA")

# Density plot----
## on log2 raw data (manually transformed)
pmexp_PC = oligo::pm(raw_data_PC); # Note that we only consider the perfect matched probes as mismatched probes are controversial
sampleNames_PC = vector()
logs_PC = vector()
for (i in 1:18)
{
  sampleNames_PC = c(sampleNames_PC,rep(as.character(pData(raw_data_PC)[i,2]),dim(pmexp_PC)[1])) #Note: here sample names in 'pData(raw_data_PC)[i,2]' must be character: 'as.character(pData(raw_data_PC)[i,2])'
  logs_PC = c(logs_PC,log2(pmexp_PC[,i]))  # log2 transformation of the raw data
}

logData_PC = data.frame(logInt=logs_PC,sampleName=sampleNames_PC)

dataHist2_PC = ggplot(logData_PC, aes(logInt, colour = sampleName)) 
dataHist2_PC + geom_density()+ ggtitle("Density plot of the log2-intensitites raw data")

## on RMA normalised data
exp_palmieri_PC <- Biobase::exprs(palmieri_eset_norm_PC)
exp_RMA_PC = exp_palmieri_PC  
sampleNames_RMA_PC = vector()
logs_RMA_PC = vector()
for (i in 1:18)
{
  sampleNames_RMA_PC = c(sampleNames_RMA_PC,rep(as.character(pData(palmieri_eset_norm_PC)[i,2]),dim(exp_RMA_PC)[1])) 
  logs_RMA_PC = c(logs_RMA_PC,exp_RMA_PC[,i])
}

logData_RMA_PC = data.frame(logInt=logs_RMA_PC,sampleName=sampleNames_RMA_PC)

dataHist2_RMA_PC = ggplot(logData_RMA_PC, aes(logInt, colour = sampleName)) 
dataHist2_RMA_PC + geom_density() + ggtitle("Density plot of the RMA calibrated prostate data")


# Box plot----

## on log2 raw data (dont need to trnasform manually as 'oligo::boxplot' performs log2 transformation by default)
oligo::boxplot(raw_data_PC, target = "core", col = FALSE,
               main = "Boxplot of log2-intensitites for Prostate raw data", xlab = "Sample Number", ylab = "Log2-intensitites", names = pData(raw_data_PC)[ , 2])

## on RMA normalised data (log2 transformed twice with the function 'oligo::boxplot', but doesnt matter I think?)
oligo::boxplot(palmieri_eset_norm_PC, target = "core", col = FALSE,
               main = "Boxplot of RMA calibrated prostate data", xlab = "Sample Number", ylab = "Log2-intensitites", names = pData(palmieri_eset_norm_PC)[ , 2]) 



##################################################################
########################### Clustering ###########################
##################################################################

# PCA----

## on log2 raw data
exp_raw_PC = log2(Biobase::exprs(raw_data_PC))
PCA_raw_PC <- prcomp(t(exp_raw_PC), scale. = FALSE) 

percentVar_PC <- round(100*PCA_raw_PC$sdev^2/sum(PCA_raw_PC$sdev^2),1)
sd_ratio_PC <- sqrt(percentVar_PC[2] / percentVar_PC[1])

dataGG_PC <- data.frame(PC1 = PCA_raw_PC$x[,1], PC2 = PCA_raw_PC$x[,2],
                        Tumor_status = pData(raw_data_PC)$Cancer.status)

ggplot(dataGG_PC, aes(PC1, PC2)) +
  geom_point(aes(colour = Tumor_status)) +
  geom_text(aes(label =pData(raw_data_PC)[ , 2] ), hjust = -0.1, vjust = -0.3) + # add point names
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar_PC[1], "%")) +
  xlim(-1000,2000)+
  ylim(-1000,1000)+ 
  ylab(paste0("PC2, VarExp: ", percentVar_PC[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio_PC) +
  scale_color_manual(values = c("#666666", "#3399FF"))


## on RMA normalised data
exp_palmieri_PC <- Biobase::exprs(palmieri_eset_norm_PC) # no need to log2 again!

PCA_RMA_PC <- prcomp(t(exp_palmieri_PC), scale = FALSE)

percentVar_RMA_PC <- round(100*PCA_RMA_PC$sdev^2/sum(PCA_RMA_PC$sdev^2),1)
sd_ratio_RMA_PC <- sqrt(percentVar_RMA_PC[2] / percentVar_RMA_PC[1])

dataGG_RMA_PC <- data.frame(PC1 = PCA_RMA_PC$x[,1], PC2 = PCA_RMA_PC$x[,2],
                            Tumor_status =  Biobase::pData(palmieri_eset_norm_PC)$Cancer.status)

ggplot(dataGG_RMA_PC, aes(PC1, PC2)) +
  geom_point(aes(colour = Tumor_status)) +
  geom_text(aes(label =pData(palmieri_eset_norm_PC)[ , 2] ), hjust = -0.1, vjust = -0.3) + # add point names
  ggtitle("PCA plot of the RMA calibrated, summarized data") +
  xlab(paste0("PC1, VarExp: ", percentVar_RMA_PC[1], "%")) +
  xlim(-200,400)+
  ylim(-200,200)+
  ylab(paste0("PC2, VarExp: ", percentVar_RMA_PC[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = sd_ratio_RMA_PC) +
  scale_color_manual(values = c("#666666", "#3399FF"))

# Dendrogram on RMA normalised data----

Biobase::pData(palmieri_eset_norm_PC)$Sample.ID_for_dendrogram = c('10A1_T','10E5_B','12A1_T','12E5_B','13A1_T','13E5_B','14A1_T','14E5_B','15A1_2_T','15D4_2_B','17A1_B','17B2_T','7A1_B','7K11_T','8A1_B','8G7_T','904_B','9A1_T')#c('10A1_T','10E5_B','12A1_T','12E5_B','13A1_T','13E5_B','14A1_T','14E5_B','15A1_2_T','15D4_2_B','17A1_B','17B2_T','7A1_B','7K11_T','8A1_B','8G7_T','904_B','9A1_T')
clustering <- exp_palmieri_PC %>% t() %>% dist(method = "euclidean") %>% 
  hclust(method = "complete")
plot(clustering, labels = (pData(palmieri_eset_norm_PC)$Sample.ID_for_dendrogram),main = "Cluster dendrogram of 18 samples")



##################################################################
########################### Annotation ###########################
##################################################################
library(affycoretools)
palmieri_eset_norm_PC <- annotateEset(palmieri_eset_norm_PC, pd.clariom.d.human) # annotating the features with annotation package; default = 'core'
head(fData(palmieri_eset_norm_PC))

length(which(!is.na(fData(palmieri_eset_norm_PC)$SYMBOL))) # there are many probes unannotated, so we chack for how many are annotated = 86161



##################################################################
#################### Intensity-based filtering ###################
##################################################################
# We filter based on means not median (???)

palmieri_means_PC <- rowMeans(Biobase::exprs(palmieri_eset_norm_PC))
hist_res_PC <- hist(palmieri_means_PC, 100, col = "#FFFFFF", freq = FALSE, 
                    main = "Histogram of the mean intensities", 
                    border = "antiquewhite4",
                    xlab = "Mean intensities") #"azure3"

## For median approach
# palmieri_means_PC_median <- rowMedians(Biobase::exprs(palmieri_eset_norm_PC))
# hist_res_PC <- hist(palmieri_means_PC_median, 100, col = "#FFFFFF", freq = FALSE, 
#                     main = "Histogram of the median intensities", 
#                     border = "antiquewhite4",
#                     xlab = "Mean intensities") #"azure3"




sh_cutoff <- shorth(palmieri_means_PC) #https://support.bioconductor.org/p/69467/#69744); shorth cutoff: 2.089

# histgram with abline
man_threshold_PC <- sh_cutoff  #!NEED DISCUSSSION WITH LORNA

hist_res_PC <- hist(palmieri_means_PC, 100, col = "#FFFFFF", freq = FALSE, 
                    main = "Histogram of the mean intensities with cutoff line",
                    border = "antiquewhite4",
                    xlab = "Mean intensities")

abline(v = man_threshold_PC, col = "coral4", lwd = 2)


filtered.set <- palmieri_eset_norm_PC[palmieri_means_PC >=sh_cutoff,] # number of probes passed the shorth cutoff
dim(filtered.set)


no_of_samples_PC <- table(paste0(pData(palmieri_eset_norm_PC)$Cancer.status)) # number of samples
no_of_samples_PC

samples_cutoff_PC <- min(no_of_samples_PC) # select the group with minimum number of samples (9 in our case) 

idx_man_threshold_PC <- apply(Biobase::exprs(palmieri_eset_norm_PC), 1,function(x){sum(x > man_threshold_PC) >=samples_cutoff_PC})
table(idx_man_threshold_PC) # number of probes that passed the short cutoff and in at least 50% samples

palmieri_manfiltered_PC <- subset(palmieri_eset_norm_PC, idx_man_threshold_PC) # filter out the probes do not meet the criteria



##################################################################
################ removing multiple mapping genes #################
##################################################################
anno_palmieri_PC <- fData(palmieri_manfiltered_PC)
anno_palmieri_PC <- subset(anno_palmieri_PC, !is.na(SYMBOL)) # removing genes with missing symbols
dim(anno_palmieri_PC)
anno_palmieri_PC[1:6,1:4]

anno_grouped_PC <- group_by(anno_palmieri_PC, PROBEID) #group by the probesetIDs 
dim(anno_grouped_PC)
anno_summarized_PC <- 
  dplyr::summarize(anno_grouped_PC, no_of_matches = n_distinct(SYMBOL)) #list out all the distinct IDs and the corresponding frequency
head(anno_summarized_PC)

anno_filtered_PC <- filter(anno_summarized_PC, no_of_matches > 1) # filter for the IDs that have mapping greater than 1
head(anno_filtered_PC)

probe_stats_PC <- anno_filtered_PC #With dim(probe_stats), we can see how many probes have been mapped to multiple genes.
nrow(probe_stats_PC)

ids_to_exlude_PC <- (featureNames(palmieri_manfiltered_PC) %in% probe_stats_PC$PROBEID) # matching the probes that need to be removed
table(ids_to_exlude_PC)

palmieri_final_PC <- subset(palmieri_manfiltered_PC, !ids_to_exlude_PC) #filter out the probes having multiple mappings; and creat the final ExpressionSet
validObject(palmieri_final_PC)

head(anno_palmieri_PC)

fData(palmieri_final_PC)$PROBEID <- rownames(fData(palmieri_final_PC))
fData(palmieri_final_PC) <- left_join(fData(palmieri_final_PC), anno_palmieri_PC)
rownames(fData(palmieri_final_PC)) <- fData(palmieri_final_PC)$PROBEID 
validObject(palmieri_final_PC)
length(which(!is.na(fData(palmieri_final_PC)$SYMBOL)))  # check for number of truly annotated probesets
boxplot(palmieri_final_PC)




##################################################################
##################  removing unannotated genes ###################
##################################################################
f_filbyNA_PC <- na.omit(fData(palmieri_final_PC)) #filter out the NA rows
ids_to_exlude_NA_PC <- (featureNames(palmieri_final_PC) %in% f_filbyNA_PC$PROBEID) #match the feature names from ExpressionSet to the annotated probesetIDs that are left; it is a logical object
palmieri_final_PC = subset(palmieri_final_PC, ids_to_exlude_NA_PC)  # subset the ExpressionSet for which includes the features that are only truly annotated
palmieri_final_PC



##################################################################
#################### removing AceView genes ######################
##################################################################

df_Aceview = fData(palmieri_final_PC)
no_of_AceView = df_Aceview[grepl("AceView", df_Aceview$GENENAME),]
head(no_of_AceView)
length(no_of_AceView$PROBEID)

f_filbyAceView_PC = df_Aceview[!grepl("AceView", df_Aceview$GENENAME),]
ids_to_match_AceView_PC <- (featureNames(palmieri_final_PC) %in% f_filbyAceView_PC$PROBEID)
palmieri_final_PC = subset(palmieri_final_PC, ids_to_match_AceView_PC)
palmieri_final_PC


##################################################################
#################### removing pseudogenes ########################
##################################################################


df_pseudo = fData(palmieri_final_PC)
no_of_pseudo = df_pseudo[grepl("pseudogene", df_pseudo$GENENAME),]
head(no_of_pseudo)
length(no_of_pseudo$PROBEID)

f_filbypseudogene_PC = df_pseudo[!grepl("pseudogene", df_pseudo$GENENAME),]
ids_to_match_pseudogene_PC <- (featureNames(palmieri_final_PC) %in% f_filbypseudogene_PC$PROBEID)
palmieri_final_PC = subset(palmieri_final_PC, ids_to_match_pseudogene_PC)
palmieri_final_PC


##################################################################
################ removing genes without GENENAME #################
##################################################################
# remove '---' genenames
df_threebar = fData(palmieri_final_PC)
no_of_threebar = df_threebar[grepl("---", df_threebar$GENENAME),]
head(no_of_threebar)
length(no_of_threebar$PROBEID)

f_filbythree_PC = df_threebar[!grepl("---", df_threebar$GENENAME),]
ids_to_match_threebar_PC <- (featureNames(palmieri_final_PC) %in% f_filbythree_PC$PROBEID)
palmieri_final_PC = subset(palmieri_final_PC, ids_to_match_threebar_PC)
palmieri_final_PC

#######################################################################################################################
################ Collapsing multiple transcript-cluster identifiers to one measure for the same genes #################
#######################################################################################################################
# Check for multiple to same gene symbol----
fData = fData(palmieri_final_PC)
n_occur <- data.frame(table(fData$SYMBOL))
n_occur<-n_occur[n_occur$Freq > 1,]
head(n_occur)

# collapse by taking gene with max mean
ids_PC = fData(palmieri_final_PC)[,c(1,3)]  ## only need the probesetIDs (1st column) and gene symbols (3rd column)
ids_PC = ids_PC[ids_PC[,2] != '',]  ## get rid of 'NA' gene symbols 86161
ID2gene = ids_PC
exprSet = exprs(palmieri_final_PC) 


{
  MAX = by(exprSet, ID2gene[,2], 
           function(x) rownames(x)[ which.max(rowMeans(x))])
  MAX = as.character(MAX)    
  exprSet = exprSet[rownames(exprSet) %in% MAX,]        
  #rownames( exprSet ) = ID2gene[ match( rownames( exprSet ), ID2gene[ , 1 ] ), 2 ] # rownames to gene symbol
}
dim(exprSet)

ids_fr_mulprocolpase_probes_PC <- (featureNames(palmieri_final_PC) %in% rownames(exprSet))
palmieri_final_PC = subset(palmieri_final_PC, ids_fr_mulprocolpase_probes_PC)
palmieri_final_PC
boxplot(palmieri_final_PC)
save(palmieri_final_PC, file = "palmieri_final_PC.RData") # This is the final ExpressionSet with 22165 expressed genes!!!!




