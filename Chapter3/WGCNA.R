# Modified from previous script 'Prostate cancer_transfer_WGCNA_v2.Rmd'

##################################################################
######################## working directory #######################
##################################################################
setwd('/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Gene identification/workspace')#Set this to your preferred working directory
##################################################################

##################################################################
################### data import for initiation ###################
##################################################################
#Pre-processed expression data----
input_palmieri_final_PC = load(file = '/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Gene identification/data_needed/palmieri_final_PC.RData')
input_palmieri_final_PC

##################################################################
######################## data preparation ########################
##################################################################
# extract and prepare the expression----
ProstateData = exprs(palmieri_final_PC) 
dim(ProstateData)
colnames(ProstateData)
boxplot(ProstateData) # sample boxplot
datExpr_PC_22165 = as.data.frame(t(ProstateData))  # transpose so row = samples; col = genes
gsg_PC = goodSamplesGenes(datExpr_PC_22165, verbose = 3); #check for genes and samples with too many missing values
gsg_PC$allOK # allgenes passed if TRUE

#Sample clustering (done this in preprocessing too)----
sampleTree_PC = hclust(dist(datExpr_PC_22165), method = "average");
plot(sampleTree_PC, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2) 

#prepare for Phenotype info (clinical traits) data ----
##Import prostate sample information
Biobase::pData(palmieri_final_PC)$Batches = as.numeric(c('1','1','1','1','2','2','2','2','2','2','2','2','1','1','1','1','1','1'))
Biobase::pData(palmieri_final_PC)$Cancer.status <- 
  c('T','B','T','B','T','B','T','B','T','B','B','T','B','T','B','T','B','T')
Biobase::pData(palmieri_final_PC)$Cancer_status = (c(1,0,1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1))
Biobase::pData(palmieri_final_PC)$CEL_ID = rownames(pData(palmieri_final_PC))
Biobase::pData(palmieri_final_PC)$Sample.ID_integrated <- 
  c('10A1_T','10E5_B','12A1_T','12E5_B','13A1_T','13E5_B','14A1_T','14E5_B','15A1_2_T','15D4_2_B','17A1_B','17B2_T','7A1_B','7K11_T','8A1_B','8G7_T','904_B','9A1_T')


traitDataPC = ((pData(palmieri_final_PC)))
dim(traitDataPC)
names(traitDataPC)
##Remove columns we dont need from the trait data frame
allTraitsPC = traitDataPC[ , -1];
dim(allTraitsPC)
names(allTraitsPC)
##Check for class
sapply(allTraitsPC, class)
sapply(allTraitsPC, is.factor)

signalSamples = rownames(datExpr_PC_22165);
traitRowsPC = match(signalSamples, allTraitsPC$CEL_ID);
datTraitsPC = as.data.frame(allTraitsPC[traitRowsPC, 5]); # use the cols that have 'numerics' for cancer status 
colnames(datTraitsPC) = "Tumour_status"
rownames(datTraitsPC) = allTraitsPC[traitRowsPC, 7];  # use the column of shorted names
collectGarbage();

table(rownames(datTraitsPC) == rownames(datExpr_PC_22165)) # check for agreement on the row names (it DOES NOT have to be agreed on the same names for further analysis)

##adding age to the datTraitsPC
datTraitsPC$age = as.numeric(c('53','53','64','64','59','59','67','67','70','70','67','67','66','66','67','67','65','65'))
##adding gleason score to the datTraitsPC
datTraitsPC$gleason_score = as.numeric(c('7','7','12','12','7','7','9','9','12','12','7','7','12','12','7','7','12','12'))


##################################################################
############## data clustering and outlier check #################
##################################################################
#Ref: http://pages.stat.wisc.edu/~yandell/statgen/ucla/WGCNA/wgcna.html
# sample network based on squared Euclidean distance
# transpose the data
A = adjacency(t(datExpr_PC_22165), type = "distance")
# this calculates the whole network connectivity
k = as.numeric(apply(A, 2, sum)) - 1
# standardized connectivity
Z.k = scale(k)

# Designate samples as outlying if their Z.k value is below the threshold
thresholdZ.k = -5  # often -2.5
# the color vector indicates outlyingness (red)
outlierColor = ifelse(Z.k < thresholdZ.k, "red", "black")
# sampleTree_flash = hclust(as.dist(1 - A), method = "average")

# Re-cluster samples
datExpr_PC_22165_rowannot <- datExpr_PC_22165
rownames(datExpr_PC_22165_rowannot) <- c('10A1_T','10E5_B','12A1_T','12E5_B','13A1_T','13E5_B','14A1_T','14E5_B','15A1_2_T','15D4_2_B','17A1_B','17B2_T','7A1_B','7K11_T','8A1_B','8G7_T','904_B','9A1_T')  # reduce the length sample name
sampleTree2_PC = hclust(dist(datExpr_PC_22165_rowannot), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors_PC = numbers2colors(datTraitsPC, signed = FALSE);

dimnames(traitColors_PC)[[2]] = paste(names(datTraitsPC))#, "C", sep = "")
datColors = data.frame(outlier = outlierColor, traitColors_PC)

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2_PC, groupLabels = names(datColors[ , 1:2]), colors = datColors[ , 1:2], # include only outlier&tumour_status
                    main = "Sample dendrogram and trait heatmap")

#Saving Expression and trait data ----
save(datExpr_PC_22165, datTraitsPC, file = "Prostate-cancer-dataInput.RData")


#####################################################################
## Step-by-step (signed) network construction and module detection ##
#####################################################################
#REF: 
##1. https://www.biostars.org/p/288153/
##2. https://support.bioconductor.org/p/73814/
##3. https://www.biostars.org/p/286227/
##4. https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/TechnicalReports/signedTOM.pdf
##5. http://pages.stat.wisc.edu/~yandell/statgen/ucla/WGCNA/wgcna.html (corrected r for Chapter12)


#Soft-thresholding power selection----
##Choose a set of soft-thresholding powers 
powers_PC = c(c(1:10), seq(from = 12, to = 20, by = 2))
##Call the network topology analysis function (IT TAKES TIME!!!)
sft_PC  = pickSoftThreshold(datExpr_PC_22165,dataIsExpr = TRUE, powerVector = powers_PC, corFnc = cor, corOptions = list(use = 'p'),networkType = "signed") #!!!for signed network

##plot results (Figure S1A)
# Create a plot
sizeGrWindow(9,5)

pdf(file = "WGCNA_sft_1.pdf", width = 13, height = 8)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft_PC$fitIndices[,1], -sign(sft_PC$fitIndices[,3])*sft_PC$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft_PC$fitIndices[,1], -sign(sft_PC$fitIndices[,3])*sft_PC$fitIndices[,2],
     labels=powers_PC,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft_PC$fitIndices[,1], sft_PC$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft_PC$fitIndices[,1], sft_PC$fitIndices[,5], labels=powers_PC, cex=cex1,col="red")
dev.off()

sft_PC$fitIndices # results for all powers 1-20
sft_PC$powerEstimate # although 'NA', but we select Power=14 from graph and 'sft_PC$fitIndices'

#Co-expression similarity and adjacency----
softPower_PC = 14;
##Checking Scale-free topology
##here we define the adjacency matrix using soft thresholding with beta=14 (IT TAKES TIME!!!)
ADJ1=abs(cor(datExpr_PC_22165,use="p"))^14
k=softConnectivity(datE=datExpr_PC_22165,power=14)
##Plot a histogram of k and a scale free topology plot
pdf(file = "WGCNA_2.pdf", width = 13, height = 8)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check Scale-free Topology\n")
dev.off() 

##adjacency matrix (IT TAKES TIME!!!)
adjacency_signed_PC = adjacency(datExpr_PC_22165, type = "signed" , power = softPower_PC); #!!!for signed network


##Topological Overlap Matrix (TOM) (IT TAKES TIME!!!)
#Turn adjacency into topological overlap
TOM_signed_PC =  TOMsimilarity(adjacency_signed_PC, TOMType="signed", verbose=5) #!!! FOR SIGNED NETWORK

dissTOM_PC = 1-TOM_signed_PC # TOM dissimilarity

#Clustering using TOM----
##Call the hierarchical clustering function
geneTree_PC = hclust(as.dist(dissTOM_PC), method = "average");
##Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree_PC, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

##We like large modules, so we set the minimum module size relatively high:
minModuleSize_PC = 50 ; #30;
##Module identification using dynamic tree cut:
dynamicMods_PC = cutreeDynamic(dendro = geneTree_PC, distM = dissTOM_PC,
                               deepSplit = 2, pamRespectsDendro = FALSE,
                               minClusterSize = minModuleSize_PC);

table(dynamicMods_PC) # No. modules and No. gene assigned

##Convert numeric labels into colors
dynamicColors_PC = labels2colors(dynamicMods_PC)
table(dynamicColors_PC)
##Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree_PC, dynamicColors_PC, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#Merging of modules whose expression profiles are very similar----
##Calculate eigengenes
MEList_PC = moduleEigengenes(datExpr_PC_22165, colors = dynamicColors_PC)
MEs_PC = MEList_PC$eigengenes
##Calculate dissimilarity of module eigengenes
MEDiss_PC = 1-cor(MEs_PC);
##Cluster module eigengenes
METree_PC = hclust(as.dist(MEDiss_PC), method = "average");
##Plot the result
sizeGrWindow(7, 6)
plot(METree_PC, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres_PC = 0.25
##Plot the cut line into the dendrogram
abline(h=MEDissThres_PC, col = "red")
##Call an automatic merging function
merge_PC = mergeCloseModules(datExpr_PC_22165, dynamicColors_PC, cutHeight = MEDissThres_PC, verbose = 3)
##The merged module colors
mergedColors_PC = merge_PC$colors;
table(mergedColors_PC)
xtable(table(mergedColors_PC)) # for latex output
##Eigengenes of the new merged modules:
mergedMEs_PC = merge_PC$newMEs;

#To see what the merging did to our module colors, we plot the gene dendrogram again, with the original and merged module colors underneath----
sizeGrWindow(12, 9)
##Plot (Figure S1B)
pdf(file = "WGCNA_3.pdf", width = 13, height = 8)

plotDendroAndColors(geneTree_PC, cbind(dynamicColors_PC, mergedColors_PC),
                    c("Dynamic Tree Cut", "Merged dynamic"), # "Merged dynamic" = "manual hybrid"
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#!In the subsequent analysis, we will use the merged module colors in mergedColors

#Saving MEs, moduleLabels, moduleColors and geneTree ----
##Rename to moduleColors
moduleColors_PC = mergedColors_PC
##Construct numerical labels corresponding to the colors
colorOrder_PC = c("grey", standardColors(50));
moduleLabels_PC = match(moduleColors_PC, colorOrder_PC)-1;
MEs_PC = mergedMEs_PC;
##Save module colors and labels for use in subsequent parts
save(MEs_PC, moduleLabels_PC, moduleColors_PC, geneTree_PC, file = "Prostate-cancer-networkConstruction-stepByStep.RData")





#####################################################################
##################### Gene network TOM visualisation ################
#####################################################################
nGenes_PC = ncol(datExpr_PC_22165);

##Gene Network Visualisation ----
##FOR aLL GENES (IT TAKES TIME!!!)
##Calculate topological overlap anew: this could be done more efficiently by saving the TOM
##calculated during module detection, but let us do it again here.
dissTOM_PC = dissTOM_PC

###FOR random 500 GENES (otherwise too big to consider 22165 genes)
nSelect_PC = 500
# For reproducibility, we set the random seed
set.seed(10);
select_PC = sample(nGenes_PC, size = nSelect_PC);
selectTOM_PC = dissTOM_PC[select_PC, select_PC];
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree_PC = hclust(as.dist(selectTOM_PC), method = "average")
selectColors_PC = moduleColors_PC[select_PC];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss_PC = selectTOM_PC^14; # power = 14
diag(plotDiss_PC) = NA;
TOMplot(plotDiss_PC, selectTree_PC, selectColors_PC, main = "Network heatmap plot, 500 genes")
# changing the heatmap color to normal
library(gplots)
pdf(file = "WGCNA_4.pdf", width = 13, height = 8)
Tomheatcol = colorpanel(250, 'red' , 'orange', 'lemonchiffon')
TOMplot(plotDiss_PC, selectTree_PC, selectColors_PC, main = "Network heatmap plot, 500 randomly selected genes", col = Tomheatcol)
dev.off()






#####################################################################
####################### Module identification #######################
#####################################################################
#Recall the objuects saved previously (these two objects were saved in both folder 'data_needed' and 'workspace', so can be recalled from either one)----
datExp_datTrait = load(file = "/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Gene identification/data_needed/Prostate-cancer-dataInput.RData")
datExp_datTrait

moduleInfo=load(file = "/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Gene identification/data_needed/Prostate-cancer-networkConstruction-stepByStep.RData")
moduleInfo

#Person correlation among the modules (by their eigengenes (i.e. first PC)) and result----
##Call the modules eigengenes 
MEs_Pcor = signif(cor(MEs_PC ,use = "p"), 2) # Person correlation of the MEs


#Diagnostic plot: module heatmap and barchart----
##module Heatmap and module barchart
##REF: https://www.jianshu.com/p/0d05d1f007e7
table(moduleColors_PC) # = table(mergedColors_PC)
unique(moduleColors_PC)
head(colnames(datExpr_PC_22165)[moduleColors_PC=="green"]) # head for probeID labelled genes in module green
featureinfo = fData(palmieri_final_PC)
head(featureinfo$SYMBOL[moduleColors_PC=="green"]) # head for gene symbol labelled genes in module green

allgenesgreen_PC = featureinfo$SYMBOL[moduleColors_PC=="green"] # list out all genes in module green
# Plot module heatmap and the eigengene all in one graph
sizeGrWindow(8,7)
which.module = "green"
ME_PC = MEs_PC[ , paste("ME", which.module, sep = "")] # here "ME" means started with letters "ME"
par(mfrow = c(2,1), mar = c(0.3, 5.5,7,2)) #c(0.3, 5.5,3,2))
plotMat(t(scale(datExpr_PC_22165[ , moduleColors_PC == which.module])),
        nrgcols = 30,
        rlabels = F,
        clabels = pData(palmieri_final_PC)$Sample.ID_integrated,
        rcols = which.module,
        #main = which.module, #先不考虑标题，更清晰一些
        cex.main = 1.2 #1
)

par(mar = c(5,4,0,0.1)) #c(5,4.2,0,0.7))
barplot(ME_PC, col = which.module, main = "", cex.main = 2, ylab = "eigengene expression", xlab = "sample" )




#Finding modules that relate to cancer status----
## Define numbers of genes and samples
nGenes_PC = ncol(datExpr_PC_22165);   
nSamples_PC = nrow(datExpr_PC_22165);  
## Recalculate module eigengenes with color labels
MEs0_PC = moduleEigengenes(datExpr_PC_22165, moduleColors_PC)$eigengenes
MEs_PC = orderMEs(MEs0_PC) # same as before 


## correlating the module eigengenes with the cancer status
ts = datTraitsPC[ , 1] # extract only tumour-status column
moduleTraitCor_PC = cor(MEs_PC, ts, use = "p"); # correlation between trait and MEs
colnames(moduleTraitCor_PC) = 'Tumour_Status'

## computing p-values for all MEs
moduleTraitPvalue_PC = corPvalueStudent(moduleTraitCor_PC, nSamples_PC); # computing p-values for all modules
colnames(moduleTraitPvalue_PC) = 'Tumour_Status'


# Module-trait heatmap (Figure S1C)----
sizeGrWindow(10,6)
## Will display correlations and their p-values
pdf(file = "WGCNA_5.pdf", width = 13, height = 8)
textMatrix_PC = paste(signif(moduleTraitCor_PC, 2), "\n(",
                      signif(moduleTraitPvalue_PC, 1), ")", sep = "");  # 'signif' for signifficance figure (s.f.)
dim(textMatrix_PC) = dim(moduleTraitCor_PC)
par(mar = c(6, 8.5, 3, 3));
## Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = as.matrix(moduleTraitCor_PC),
               xLabels = 'Tumour_Status',#names(datTraitsPC),
               yLabels = names(MEs_PC),
               ySymbols = names(MEs_PC),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),     #colors = blueWhiteRed(50)
               textMatrix = textMatrix_PC,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off() 




# Barchart of module significance as average (mean) gene significance (Figure S1D)----
Cancer_stat = as.data.frame(datTraitsPC$Tumour_status);
names(Cancer_stat) = "Tumour_status"

# Computing GS
GS1 = as.numeric(cor(Cancer_stat, datExpr_PC_22165, use = "p")); # same values as 'geneTraitSignificance_PC'; BUT one 'numerics' one 'data.frame'; GS is the correlation r value!

# taking absolute values of GS
GeneSignificance_abs = abs(GS1)

# taking average
ModuleSig = tapply(GeneSignificance_abs , moduleColors_PC, mean, na.rm = T)

# plot Module significance bar plot
sizeGrWindow(8,7)

pdf(file = "WGCNA_7.pdf", width = 13, height = 8)
par(mfrow = c(1,1))

plotModuleSignificance(GeneSignificance_abs, moduleColors_PC)
dev.off()


#Eigengenes Network Visualisation----
# Eigengenes Network Visualisation
# Recalculate module eigengenes
MEs_PC = moduleEigengenes(datExpr_PC_22165, moduleColors_PC)$eigengenes
# Isolate Cancer.status from the clinical traits
Cancer_stat = as.data.frame(datTraitsPC$Tumour_status);
names(Cancer_stat) = "tumor"
# Add the weight to existing module eigengenes
MET_PC = orderMEs(cbind(MEs_PC, Cancer_stat))

# Plot the relationships among the eigengenes and the trait
par(cex = 0.9)
plotEigengeneNetworks(MET_PC, "", marDendro = c(1,9,5,4), marHeatmap = c(6,8,1,0.8), cex.lab = 0.8, xLabelsAngle
                      = 90)




#####################################################################
#### Gene relationship to trait and important modules (GS vs MM) ####
#####################################################################
#Gene relationship to trait and important modules
#Here, we calculate the **Module Membership (MM)** and **Gene Significance (GS)**

#Gene relationship to trait and important modules: GS and MM ----
Cancer_stat = as.data.frame(datTraitsPC$Tumour_status);
names(Cancer_stat) = "Tumour_status"
moduleTraitCor_PC[ , "Tumour_Status"] # look for correlation of the module eigengene and clinical trait


#Module Membership----
##First, calculate the correlation matrix between genes and corresponding modules
##names (colors) of the modules
modNames_PC = substring(names(MEs_PC), 3)
geneModuleMembership_PC = as.data.frame(cor(datExpr_PC_22165, MEs_PC, use = "p"));

# Calculate the Person correlation coefficient matrix of each module and its genes
MMPvalue_PC = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_PC), nSamples_PC));

names(geneModuleMembership_PC) = paste("MM", modNames_PC, sep="");
names(MMPvalue_PC) = paste("p.MM", modNames_PC, sep="");


#Gene Significance, the same values as GS1 calculated for mod sig barplot----
# Second, calculate the correlation matrix between trait and genes; only continuous trait can be computed. If the trait is discrete, need to convert first to 1/0 (1=cancer; 0=normal)
geneTraitSignificance_PC = as.data.frame(cor(datExpr_PC_22165, Cancer_stat, use = "p"));
GSPvalue_PC = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_PC), nSamples_PC));

names(geneTraitSignificance_PC) = paste("GS.", names(Cancer_stat), sep="");
names(GSPvalue_PC) = paste("p.GS.", names(Cancer_stat), sep="");


#Intramodular analysis: Identifying genes with high GS and MM ---- 

# For GREEN module
module = "green"
column = match(module, modNames_PC);
moduleGenes_PC = moduleColors_PC==module;
# Open a pdf file

# Create a plot (Figure S1E)
sizeGrWindow(7, 7);
pdf(file = "WGCNA_8.pdf", width = 13, height = 8)
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership_PC[moduleGenes_PC, column]),
                   abs(geneTraitSignificance_PC[moduleGenes_PC, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for cancer status(malignant and normal)",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)+ abline(h = 0.3, col = 'red', lwd =2) + abline(v = 0.9, col = 'red', lwd =2)

dev.off() 

#####################################################################
########### Summary output of network analysis results ##############
#####################################################################

# Output of probe IDs of all genes
names(datExpr_PC_22165)
# Output of probe IDs of the genes belong to the orangered4 or turquoise modules
names(datExpr_PC_22165)[moduleColors_PC=="green"]

#Change the following to your working directory
annot_PC <- read.csv("/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Gene identification/data_needed/GeneAnnotation_for_WGCNA.csv")
dim(annot_PC)
names(annot_PC)
probes_PC = names(datExpr_PC_22165)
probes2annot_PC = match(probes_PC, annot_PC$PROBEID)
# The following is the number or probes without annotation:
sum(is.na(probes2annot_PC))
# Should return 0

# Now we create a new data frame containing the following information for all of the probes: 
# Probe ID, Public Gene ID, Gene Symbol, Description, Group, Module Colour, Gene Significance for T2D, Module Membership and P-values in all modules. 
# The data file ‘geneInfo_Islet_power10.csv’ was exported.

# Create the starting data frame
geneInfo0_PC = data.frame(ProbeIDs = probes_PC,
                          Gene_Symbol = annot_PC$SYMBOL[probes2annot_PC],
                          EntrezID = annot_PC$ENTREZID_wgcna[probes2annot_PC],
                          
                          moduleColor = moduleColors_PC,
                          geneTraitSignificance_PC,
                          GSPvalue_PC)                              
# Order modules by their significance for disease status
modOrder_PC = order(-abs(cor(MEs_PC, Cancer_stat, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership_PC))
{
  oldNames_PC = names(geneInfo0_PC)
  geneInfo0_PC = data.frame(geneInfo0_PC, geneModuleMembership_PC[, modOrder_PC[mod]],
                            MMPvalue_PC[, modOrder_PC[mod]]);
  names(geneInfo0_PC) = c(oldNames_PC, paste("MM.", modNames_PC[modOrder_PC[mod]], sep=""),
                          paste("p.MM.", modNames_PC[modOrder_PC[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder_PC = order(geneInfo0_PC$moduleColor, -abs(geneInfo0_PC$GS.Tumour_status));
geneInfo_PC = geneInfo0_PC[geneOrder_PC, ]
write.csv(geneInfo_PC, file = "moduleInfo_allgenes_prostate_power14.csv")  # summarization for all genes 22165


# For genes in GREEN module only
gene_green_PC= geneInfo_PC[geneInfo_PC$moduleColor == 'green','ProbeIDs']
ids_green_PC = geneInfo_PC$ProbeIDs %in% gene_green_PC
gene_green_list_PC= subset(geneInfo_PC, ids_green_PC)
write.csv(gene_green_list_PC,"moduleInfo_green_prostate_power14.csv",quote = F)





#####################################################################
####################### hub gene selection ##########################
#####################################################################
#REF:
#https://mp.weixin.qq.com/s/VKOc5LpcWDuG_wZgBO1s2w
#https://www.jianshu.com/p/0d05d1f007e7

# Recalculate topological overlap if needed
TOM_signed_PC = TOM_signed_PC;
# Read in the annotation file
annot_PC <- read.csv("/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Gene identification/data_needed/GeneAnnotation_for_WGCNA.csv")

# Select modules
modules_PC = "green" #c("brown", "red");
# Select module probes
probes_PC = colnames(datExpr_PC_22165)
# e =(moduleColors_PC == modules_PC) #!!!the SAME as 'is.finite(match(moduleColors_PC, modules_PC));' below
# y = probes_PC[e]
# head(y)
inModule_PC = is.finite(match(moduleColors_PC, modules_PC));
modProbes_PC = probes_PC[inModule_PC];  # !!! same probeID as in the 'gene_green_list_PC', BUT not in the same order
head(modProbes_PC)
#https://www.biostars.org/p/335283/ :
# top.genes.logical = cbind(modProbes_PC,top) ;
# top.genes.sort = as.matrix (top.genes.logical[order(-top, modProbes_PC)])
# modGenes_PC = annot_PC$SYMBOL[match(top.genes.sort, annot_PC$PROBEID)];#annot_PC$SYMBOL[match(modProbes_PC, annot_PC$PROBEID)];
ids_greenProbe_PC <- (probes_PC %in% modProbes_PC) #match the feature names from ExpressionSet to the annotated probesetIDs that are left; it is a logical object
greenMod_annot_PC = subset(fData(palmieri_final_PC), ids_greenProbe_PC)  # subset the ExpressionSet for which includes the features that are only truly annotated
#View(greenMod_annot_PC) #Now, same order as 'modProbes_PC'
head(greenMod_annot_PC) #Now, same order as 'modProbes_PC'; We can use this data frame for CHECKING if the top ranking IMC genes are haveing the correct 'gene symbol'!
modGenes_PC = annot_PC$SYMBOL[match(modProbes_PC, annot_PC$PROBEID)];# the same as 'greenMod_annot_PC' !!!
# Select the corresponding Topological Overlap
modTOM_PC = TOM_signed_PC[inModule_PC, inModule_PC];
dimnames(modTOM_PC) = list(modProbes_PC, modProbes_PC)


# REF:
#https://www.jianshu.com/p/0d05d1f007e7

#(1)calculating the Intramodular connectivity ----
##! For signed network
# b2 = intramodularConnectivity.fromExpr(datExpr_PC_22165, moduleColors_PC,networkType = "signed", power = 14)
# head(b2)
Alldegrees1_signed_PC = intramodularConnectivity(adjacency_signed_PC,moduleColors_PC )
head(Alldegrees1_signed_PC)
#b2 = Alldegrees1_signed_PC

#(2)calculating relationship between GS and intramodular connectivity ----
module = "green"
P = as.data.frame(datTraitsPC$Tumour_status) # specific clinical info
names(P) = "Tumour_Status"
GS1 = as.numeric(cor(P, datExpr_PC_22165, use = "p"))
GeneSignificance = abs(GS1)

# Single module GS vs IMC plot only (Figure S1F)----
pdf(file = "WGCNA_9.pdf", width = 13, height = 8)
which.color="green"; 
restrictGenes=moduleColors_PC==which.color 
verboseScatterplot(Alldegrees1_signed_PC$kWithin[ restrictGenes],
                   GeneSignificance[restrictGenes],
                   col=which.color,
                   xlab="Intramodular Connectivity",
                   ylab="Gene Significance",
                   abline = TRUE,
                   main = 'Gene Significance vs. Intramodular Connectivity (Green Module)')
dev.off()


#(3) Calculating the connectivity of all genes in the module. Screen the hub genes----
#!!! Calculating module membership(MM); or generalising IMC for all genes on the array
datKME_PC = signedKME(datExpr_PC_22165, MEs_PC, outputColumnName = "MM.")
head(datKME_PC)  

write.csv(datKME_PC, file = "datKME_all_PC.csv")  # gene connectivity in the module



# plot the MM vs IMC separately 
pdf(file = "WCGNA_10.pdf", width = 13, height = 8)
#power of MM^5
which.color="green"; 
restrictGenes=moduleColors_PC==which.color 
verboseScatterplot(Alldegrees1_signed_PC$kWithin[ restrictGenes],
                   abs(datKME_PC[restrictGenes, paste("MM.", which.color, sep="")])^5, col=which.color,
                   xlab="Intramodular Connectivity",
                   ylab="abs(Module Membership)^5",
                   abline = TRUE,
                   main = 'Module Membership vs. Intramodular Connectivity (Green Module)')
dev.off()

# (Figure S1G)
pdf(file = "WGCNA_11.pdf", width = 13, height = 8)
#power of MM^1
which.color="green";  
restrictGenes=moduleColors_PC==which.color 
verboseScatterplot(Alldegrees1_signed_PC$kWithin[ restrictGenes],
                   abs(datKME_PC[restrictGenes, paste("MM.", which.color, sep="")])^1, col=which.color,
                   xlab="Intramodular Connectivity",
                   ylab="abs(Module Membership)",
                   abline = TRUE,
                   main = 'Module Membership vs. Intramodular Connectivity (Green Module)')
dev.off()


# export the results from 'softConnectivity()' and 'intramodularConnectivity()' so we can checking the ranking of the top genes
write.csv(IMConn, file = "IMConn_by_softCon.csv") 
write.csv(Alldegrees1_signed_PC$kWithin[ restrictGenes], file = "IMC_by_intramod.csv") 


# TOP5% topgenes:  KME  measure ====
TOM_IMC_green <- subset(Alldegrees1_signed_PC,rownames(Alldegrees1_signed_PC) %in% greenMod_annot_PC$PROBEID)
TOM_IMC_green <- TOM_IMC_green[order(-TOM_IMC_green$kWithin) , ]  # 'order(-)' to order from large to small NOT '-order()'
topgenes_TOM_IMC <- TOM_IMC_green[1:220 , ]

write.csv(topgenes_TOM_IMC,"topgenes_TOM_IMC_220.csv")  #！！！no need to convert to gene symbol， use probeid and input to Venndiagram software

# TOP5% topgenes:  tom_imc measure ====
MM_IMC_green <- subset(datKME_PC,rownames(datKME_PC) %in% greenMod_annot_PC$PROBEID)
MM_IMC_green <- MM_IMC_green[order(-abs(MM_IMC_green$MM.green)) , ]  # 'order(-)' to order from large to small NOT '-order()'
topgenes_KME_IMC <- MM_IMC_green[1:220 , ]
min(MM_IMC_green$MM.green)
min(topgenes_KME_IMC$MM.green)

write.csv(topgenes_KME_IMC,"topgenes_KME_IMC_220.csv")


topgenes_overlap <- subset(topgenes_KME_IMC ,rownames(topgenes_KME_IMC) %in% rownames(topgenes_TOM_IMC)) # 165 overlapped genes
min(topgenes_overlap$MM.green) #minimum MM among the 165 hub genes

# gene symbols
topgenes_overlap_KME_tomIMC <- subset(fData(palmieri_final_PC), rownames(fData(palmieri_final_PC)) %in% rownames(topgenes_overlap))
topgenes_overlap_KME_tomIMC <- topgenes_overlap_KME_tomIMC$SYMBOL

write.csv(topgenes_overlap_KME_tomIMC,"topgenes_overlap_KME_tomIMC.csv")



# GS check for the overlap genes
topgenes_overlap_for_GScheck <- subset(fData(palmieri_final_PC), rownames(fData(palmieri_final_PC)) %in% rownames(topgenes_overlap))

topgenes_overlap_GS_1 <- subset(geneTraitSignificance_PC, rownames(geneTraitSignificance_PC) %in% topgenes_overlap_for_GScheck$PROBEID) # 'geneTraitSignificance_PC' has '-' values as NOT taking absolute values
topgenes_overlap_GS_1$PROBEID <- rownames(topgenes_overlap_GS_1)

topgenes_overlap_GS_2 <- subset(topgenes_overlap_for_GScheck, topgenes_overlap_for_GScheck$PROBEID %in% rownames(geneTraitSignificance_PC))

topgenes_overlap_GS <- merge(topgenes_overlap_GS_1,topgenes_overlap_GS_2 ,by.y='PROBEID',by.x='PROBEID')
min(topgenes_overlap_GS$GS.Tumour_status) # minimum GS among the 165 hub genes























