# Working directory----
setwd("/Users/zhuofanmou/Documents/PhD Exeter/ClariomD/alternative_splicing/R_workspace_v2/limmaVoom_DGEA_TCGA")
getwd()



# Ref:
#https://support.bioconductor.org/p/77664/
#LIMMA user guide



# Library ----

library(limma)
library(edgeR)


# read in raw counts data----
x <- read.table("/Users/zhuofanmou/Documents/PhD Exeter/ClariomD/alternative_splicing/R_workspace_v2/limmaVoom_DGEA_TCGA/datExp.counts.prad.txt",header = T,sep="\t",row.names = 1,check.names = F) # must add check.names = F, otherwise the '-' wont display for sample name 'xxxx-xxxx-xxxx'
dim(x) # 56602*551
View(x)



# pheno Info
pheno = read.table("target.txt",header = T,sep="\t",row.names = 1)
View(pheno)

pheno <- pheno[match(colnames(x), rownames(pheno)),] 
table(colnames(x) == rownames(pheno))

#Design matrix (NOTE for LIMMAVOOM a design matrix is needed!!!)
Group = factor(pheno$Group,levels=c('Tumour','Normal'))
design = model.matrix(~0+Group)
colnames(design) <- c('Tumour','Normal')

# design = model.matrix(~0+Group+pheno$Sex+pheno$Age)
# colnames(design) <- c('OA','normal','Sex','Age')
design


# Gordon suggested to use either TMM (in the calcNormFactors) or not to use TMM but do the quantile normalisation (in voom)
# voom by default do not do the QN (i.e., normalise.method = "none" or dont put anything); so if we want to do QN, we need to set in voom, BUT dont use TMM in the calcNormFactors anymore (i.e., method = "none")!!!

# Normalisation and filtering 

y <- DGEList(counts = x)
# View(y$samples)
# y$samples$group <- pheno$Group # the feeding of group is crucial and be careful that we must use y$counts in the 'voom' function!

# boxplot((y$counts),range=0,ylab="log2 intensity")



# # option1: to filter genes before normalization! BASED on LIMMA userguide!
# A <- rowSums(y$counts)
# isexpr <- A > 50
# 
# y <- y[isexpr , , keep.lib.size = FALSE]
# dim(y)



# # option2: to filter genes before normalization! BASED on LIMMA userguide!
#isexpr <- rowSums(cpm(y, log = TRUE) > 3) >= 1# according to the paper for GSe114007
# isexpr <- rowSums(cpm(y) > 1) >= 3
# y <- y[isexpr,,keep.lib.sizes = FALSE]
# dim(y)

# # option3: to filter genes before normalization! BASED on LIMMA userguide!

keep <- filterByExpr(y, design)  # Note!!! in order to do this, we must pre-define the 'group' in 'y$samples'
y <- y[keep ,, keep.lib.sizes = FALSE]
dim(y)  #24957*551

#It is usual to apply scale normalization to RNA-seq read counts, and the TMM normalization method [33] in particular has been found to perform well in comparative studies. This can be applied to the DGEList object:
# Normalisation with TMM (this following code still then need to be feeded into voom BUT not for normalized.method ="quantile")
y <- calcNormFactors(y, method = "TMM")
#y <- calcNormFactors(y, method = "none")

View(y$counts)

# Normalisation with quantile using VOOM dunction (using the counts matrix y$counts as input NOT the (TMM) DElist y)
#v <- voom(y$counts, normalize.method = "quantile")  # if the sample group was defined (for option3), here, we must use the count matrix (y$counts) instead of DElist y
v <- voom(y,design,plot=TRUE)

View(v$E)
dim(v)

boxplot((v$E),range=0,ylab="log2 intensity")

# HERE we switch v back to y just for the sake of folloing scripts and keep the workflow consistent

y <- v
class(y)
View(y$E)

# output data
exp_norm_wtFilter <- y$E
write.table(exp_norm_wtFilter,"datExp.normalized.txt",row.names = T,col.names = T,quote = F,sep="\t")
save(exp_norm_wtFilter, file = "datExp.normalized.rda") # save the unfiltered normalised expression matrix










# Annotating in y
anno = read.csv("Annotation.csv",header = T) # Gene symbol is not label as 'NA' but "", so we need 'na.strings = c("","NA")'

anno<- sapply(anno, as.character)
anno[is.na(anno)] <- ""
anno <- as.data.frame(anno)
class(anno)
str(anno)
View((anno))
colnames(anno) = c("ID","Gene.Symbol")
View(head(anno))


annotLookup <- anno[which(anno$ID %in% rownames(y)),]
table(annotLookup$ID==rownames(y))
View(annotLookup)
annotLookup <- annotLookup[match(rownames(y), annotLookup$ID),] # match annotLookup ID (thus whole annotLookup data frame) to the same order as rownames in y
table(rownames(y) == annotLookup$ID) # check that annots are aligned
y$ID <- annotLookup$ID
# y$EntrezID <- annotLookup$EntrezID
y$Gene.Symbol <- annotLookup$Gene.Symbol
View(y$E)
View(y$Gene.Symbol)
# # y$genes$ensembl_gene_id <- annotLookup$ensembl_gene_id
# # y$genes$entrezgene <- annotLookup$entrezgene
# # y$genes$gene_biotype <- annotLookup$gene_biotype
# # y$genes$external_gene_name <- annotLookup$external_gene_name
# 
# # Probes filtering
# Control <- y$genes$ControlType == 1L; table(Control) # remove control probes/reference probes
# #NoSymbol <- is.na(y$genes$Gene.Symbol);table(NoSymbol)  # remove unannotated probes
# #NoSymbol <- (y$genes$Gene.Symbol == "");table(NoSymbol)  # remove unannotated probes; optional!!!
# IsExpr <- rowSums(y$other$rIsWellAboveBG > 0) >= 12;table(IsExpr) # remove non-expressed probes: row sums of well above BG less than 50% of the samples
# #Isdup <- duplicated(y$genes$GeneName);table(Isdup) # Optional!!!
# 
# # yfilt <- y[!Control & !NoSymbol & IsExpr & !Isdup, ] 
# yfilt <- y[!Control  & IsExpr, ] 
# 
# 
# dim(y) # 62976*22
# dim(yfilt) # 50771*22

# check for duplicates
table(duplicated(y$ID)) # indicating no one-to-many case
# table(duplicated(y$EntrezID))
table(duplicated(y$Gene.Symbol))

# # for replicate probes, replace values with the mean
# # ID is used to identify the replicates
# yfilt <- avereps(yfilt, ID =yfilt$genes$AgilentID )
# dim(yfilt) # 48305*22

y$Gene.Symbol[is.na(y$Gene.Symbol)] <- ""


# collapse probes of many-to-one gene symbol by taking gene with max mean----
ids = annotLookup [ , c(1,2)]  ## only need the probesetIDs (1st column) and gene symbols (3rd column)
View(ids)
table(is.na(ids$Gene.Symbol)) # it has 'NA's!!!
table(ids$Gene.Symbol == "")  # No ''
#table(ids$Gene.Symbol == "---")  # No '---'
ids = ids[ids[,2] != '',]  ## get rid of 'NA' gene symbols 86161
ids = na.omit(ids)
ID2gene = ids

y<- y[match( ids$ID , rownames(y)),] # match annotLookup ID (thus whole annotLookup data frame) to the same order as rownames in y
dim(y)
exprSet = y$E
table(ID2gene$ID == rownames(exprSet))

{
  MAX = by(exprSet, ID2gene[,2],
           function(x) rownames(x)[ which.max(rowMeans(x))])
  MAX = as.character(MAX)
  exprSet = exprSet[rownames(exprSet) %in% MAX,]
  #rownames( exprSet ) = ID2gene[ match( rownames( exprSet ), ID2gene[ , 1 ] ), 2 ] # rownames to gene symbol
}

dim(exprSet)
View(exprSet)

y<- y[match( rownames(exprSet) , rownames(y)),] # match annotLookup ID (thus whole annotLookup data frame) to the same order as rownames in y


dim(y)
View(y)
table(rownames(exprSet) == rownames(y))

# re-annotate it again!
annotLookup <- anno[which(anno$ID %in% rownames(y)),]
table(annotLookup$ID==rownames(y))
View(annotLookup)
annotLookup <- annotLookup[match(rownames(y), annotLookup$ID),] # match annotLookup ID (thus whole annotLookup data frame) to the same order as rownames in y
table(rownames(y) == annotLookup$ID) # check that annots are aligned
y$ID <- annotLookup$ID
#y$EntrezID <- annotLookup$EntrezID
y$Gene.Symbol <- annotLookup$Gene.Symbol
View(y)
##########################################################################################

#Check for duplicates
table(duplicated(y$ID)) # all unique now
table(duplicated(y$Gene.Symbol)) # still duplicates gene symbol, indicating many-to-one case
NoSymbol <- (y$Gene.Symbol == "");table(NoSymbol)  # No. unannotated probes



exp_filtered = y$E
length(y$ID)
length(y$Gene.Symbol)
# rownames(exp_filtered) <- y$ID
# colnames(exp_filtered)= str_extract(colnames(exp_filtered),"GSM\\d*")
table(duplicated(rownames(exp_filtered)))

exp_filtered[1:2,1:2]
View(exp_filtered)

boxplot(exp_filtered)


pdf(file = "BoxplotAfterNorm_filtered.pdf", width = 20, height = 8)

# boxplot
par(mar=c(6,3,2,1))
boxplot(exp_filtered,las=2,pch=16,
        cex=0.5,cex.axis=0.8,
        notch=T,outline = F,
        main=("Boxplot of the TCGA-PRAD"))

# draw median value levels
abline(h = mean(apply(exp_filtered,2,median)), col = "red", lwd=2)  # median
abline(h = mean(apply(exp_filtered,2,quantile)[2,]), col = rgb(0,0,255,80,maxColorValue = 255), lwd=2, ) # lower quantile
abline(h = mean(apply(exp_filtered,2,quantile)[4,]), col = rgb(0,0,255,80,maxColorValue = 255), lwd=2)  # upper quantile

dev.off()



# output data
write.table(exp_filtered,"datExp.normalized_filtered.txt",row.names = T,col.names = T,quote = F,sep="\t")
save(exp_filtered, file = "datExp.normalized_filtered.rda") # same the unfiltered normalised expression matrix
























