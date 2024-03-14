

##################################################################################################################################################################################
################################################################################## Set-up ########################################################################################
##################################################################################################################################################################################

# Working directory----
setwd("/Users/zhuofanmou/Documents/Phd/papers/paper3/data/GeneMeta_prostate_tissue/meta_analysis/differential_v2")
getwd()


# Libraries
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")



library(limma)
library(ggplot2)
library(Biobase)
library(arrayQualityMetrics)

options(stringsAsFactors = F)
library(DLL)
library(stringr)
library(ggsci)
library(ggplot2)
library(Rtsne)
library(DLL)
library(reticulate)
library(umap)
library(cowplot)

options(stringsAsFactors = F)
library(stringr)
library(limma)
library(cowplot)
library(openxlsx)
library(ggpubr)
library(ggthemes)
library(dplyr)


#-----------------------------------------------------------------##
#              DEG analysis per dataset                           ##
#-----------------------------------------------------------------##   

##################################################################################################################################################################################
############################################################################ TCGA ################################################################################################
##################################################################################################################################################################################


# read in ExpressionSet

TCGA <- readRDS('TCGA-PRAD_eSet.RDS')
TCGA
View(exprs(TCGA))
View(pData(TCGA))
View(fData(TCGA))
boxplot(TCGA)


# Expression
exp_TCGA <- exprs(TCGA)
View(exp_TCGA)
dim(exp_TCGA)
write.table(exp_TCGA, file="exp_TCGA.txt", sep="\t", quote=F, row.names=F)

# Clinical
clinic_TCGA <- pData(TCGA)
dim(clinic_TCGA)
write.table(clinic_TCGA, file="target_TCGA.txt", sep="\t", quote=F, row.names=F)


# Annotation
anno <- readRDS("PCaDB_Gene_Annotation.RDS")
anno_spec <- anno[ anno$gene_name %in% c("ZWINT", 
                                         "APEX1",
                                         "APOF",
                                         "BZW2",
                                         "CACNA1D",
                                         "CAMKK2",
                                         "DANCR",
                                         "DDT",
                                         "DDTL",
                                         "DPM2",
                                         "EEF1D",
                                         "ERBB3",
                                         "FASN",
                                         "FBP1",
                                         "FKBP4",
                                         "GJB1",
                                         "GLT8D2",
                                         "GMDS",
                                         "HOXC6",
                                         "IMPDH2",
                                         "LAMTOR2",
                                         "MRPL24",
                                         "MYC",
                                         "MYL5",
                                         "NOP16",
                                         "OR51F2",
                                         "PGM5P3-AS1",
                                         "PHB2",
                                         "PILRB",
                                         "POLD2",
                                         "PSMG4",
                                         "PTPRF",
                                         "RAB25",
                                         "RPLP0",
                                         "SDF4",
                                         "SLC7A11",
                                         "SMIM22",
                                         "SMUG1",
                                         "SNHG19",
                                         "SNORD49A",
                                         "TMED3",
                                         "TMEM183A",
                                         "TMEM97",
                                         "TPI1",
                                         "TRPM4",
                                         "TRPV6",
                                         "TTLL12",
                                         "ZNHIT1"), ] # for genes of interest only
View(anno_spec)
dim(anno_spec) #48*27
  


# Processing for single genes
pheno_TCGA = read.table("target_TCGA.txt",header = T,sep="\t")
rownames(pheno_TCGA) <- pheno_TCGA$sample_id
View(pheno_TCGA)  
unique(pheno_TCGA$sample_type)
table(pheno_TCGA$sample_type)
table(colnames(TCGA) == rownames(pheno_TCGA)) # 547 samples



# Assigning gene names to the rows of expression matrix
same=intersect(row.names(anno_spec),row.names(exp_TCGA)) #以'df_clinic.genemarkers' 的行名来找因为它在前而且他是最一开始的行名顺序，与‘combind_model...’一致的

genes_anno=anno_spec[same,] # annotation for selected genes
View(genes_anno)
dim(genes_anno)
head(genes_anno)
class(genes_anno)
write.table(genes_anno, file="genes_anno_TCGA.txt", sep="\t", quote=F, row.names=F)

genes_exp <- exp_TCGA[same,] # expression for selected genes
genes_exp <- as.data.frame(genes_exp)

##Check if row names matched
table(row.names(genes_exp) == row.names(genes_anno))



# Convert row names into gene names
genes_exp$geneName <- genes_anno$gene_name
#View(as.data.frame(genes_exp[ , 548]))
rownames(genes_exp) <- genes_exp$geneName
genes_exp <- genes_exp[ , -ncol(genes_exp)]
View(genes_exp)
dim(genes_exp)
head(genes_exp)
class(genes_exp)
##


# Transpose Expression matrix
data <- t(genes_exp)
View(data)
class(data)


# Extract clinical information
class(pheno_TCGA)
Type <- pheno_TCGA[ , colnames(pheno_TCGA) %in% c("sample_id", "sample_type")] 
View(Type)

same_type=intersect(row.names(data),row.names(Type)) #以'df_clinic.genemarkers' 的行名来找因为它在前而且他是最一开始的行名顺序，与‘combind_model...’一致的
length(same_type)
df1=data[same_type,] # expression for selected genes
View(df1)
dim(df1)
head(df1)
class(df1)

df2 <- Type[same_type,] # pheno of the samples
View(df2)
class(df2)

##Check if row names matched
table(row.names(df1) == row.names(df2))

exp_TCGA_genes <- cbind(df1,df2)
View(exp_TCGA_genes)


#输出目标基因的表达量
outTab=exp_TCGA_genes
names(outTab)[ncol(outTab)] <- "Type"
View(outTab)
outTab=cbind(ID=row.names(outTab), outTab)
write.table(outTab, file="geneExp_TCGA.txt", sep="\t", quote=F, row.names=F)

#Differential analysis

# LOOP for diff analysis

genes <- rownames(genes_exp)

for(gene in genes) {

exp <- outTab [ , c(gene,"Type")] 
View(exp)
str(exp)
group=levels(factor(exp$Type))
exp$Type=factor(exp$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#create boxplot
boxplot=ggboxplot(exp, x="Type", y=gene, color="Type",
                  xlab="",
                  ylab=paste0(gene, " expression"),
                  legend.title="Type",
                  palette = c("blue","red"),
                  add = "jitter")+ 
  stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")

#save the plot
pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
print(boxplot)
dev.off()

}


# ## "ZWINT"
# gene <- "ZWINT"
# exp <- outTab [ , c("ZWINT","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="ZWINT", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ## "PCBP1"
# 
# gene <- "PCBP1"
# exp <- outTab [ , c("PCBP1","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="PCBP1", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# ## "PABPN1" 
# 
# gene <- "PABPN1"
# exp <- outTab [ , c("PABPN1","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="PABPN1", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ## "PTPRF"
# gene <- "PTPRF"
# exp <- outTab [ , c("PTPRF","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="PTPRF", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# 
# ## "DANCR"
# 
# gene <- "DANCR"
# exp <- outTab [ , c("DANCR","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="DANCR", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# ## "MYC"
# 
# 
# gene <- "MYC"
# exp <- outTab [ , c("MYC","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="MYC", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ## "EEF1D"
# 
# 
# gene <- "EEF1D"
# exp <- outTab [ , c("EEF1D","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="EEF1D", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# ## "FBP1"
# 
# gene <- "FBP1"
# exp <- outTab [ , c("FBP1","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="FBP1", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ##"LAMTOR2"
# 
# 
# 
# gene <- "LAMTOR2"
# exp <- outTab [ , c("LAMTOR2","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="LAMTOR2", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()





##################################################################################################################################################################################
################################################################### GSE21034 (Taylor/MSKCC) ######################################################################################
##################################################################################################################################################################################

GSE21034 <- readRDS('Taylor_eSet.RDS')
GSE21034
View(exprs(GSE21034))
View(pData(GSE21034))
View(fData(GSE21034))
boxplot(GSE21034)


# Expression
exp_GSE21034 <- exprs(GSE21034)
View(exp_GSE21034)
dim(exp_GSE21034)
write.table(exp_GSE21034, file="exp_GSE21034.txt", sep="\t", quote=F, row.names=F)

# Clinical
clinic_GSE21034 <- pData(GSE21034)
dim(clinic_GSE21034)
write.table(clinic_GSE21034, file="target_GSE21034.txt", sep="\t", quote=F, row.names=F)


# Annotation
anno <- readRDS("PCaDB_Gene_Annotation.RDS")
anno_spec <- anno[ anno$gene_name %in% c("ZWINT", 
                                         "APEX1",
                                         "APOF",
                                         "BZW2",
                                         "CACNA1D",
                                         "CAMKK2",
                                         "DANCR",
                                         "DDT",
                                         "DDTL",
                                         "DPM2",
                                         "EEF1D",
                                         "ERBB3",
                                         "FASN",
                                         "FBP1",
                                         "FKBP4",
                                         "GJB1",
                                         "GLT8D2",
                                         "GMDS",
                                         "HOXC6",
                                         "IMPDH2",
                                         "LAMTOR2",
                                         "MRPL24",
                                         "MYC",
                                         "MYL5",
                                         "NOP16",
                                         "OR51F2",
                                         "PGM5P3-AS1",
                                         "PHB2",
                                         "PILRB",
                                         "POLD2",
                                         "PSMG4",
                                         "PTPRF",
                                         "RAB25",
                                         "RPLP0",
                                         "SDF4",
                                         "SLC7A11",
                                         "SMIM22",
                                         "SMUG1",
                                         "SNHG19",
                                         "SNORD49A",
                                         "TMED3",
                                         "TMEM183A",
                                         "TMEM97",
                                         "TPI1",
                                         "TRPM4",
                                         "TRPV6",
                                         "TTLL12",
                                         "ZNHIT1"), ] # for genes of interest only






View(anno_spec)
dim(anno_spec) #48*27



# Processing for single genes
pheno_GSE21034 = read.table("target_GSE21034.txt",header = T,sep="\t")
rownames(pheno_GSE21034) <- pheno_GSE21034$sample_id
View(pheno_GSE21034)  
unique(pheno_GSE21034$pcadb_group)
table(pheno_GSE21034$pcadb_group)
table(colnames(GSE21034) == rownames(pheno_GSE21034)) # Note: the PCaDB annotated the sample id differently, wheras the original expression data has the sample id with GSMxxx from the GEO



# Assigning gene names to the rows of expression matrix
same=intersect(row.names(anno_spec),row.names(exp_GSE21034)) #以'df_clinic.genemarkers' 的行名来找因为它在前而且他是最一开始的行名顺序，与‘combind_model...’一致的

genes_anno=anno_spec[same,] # annotation for selected genes
View(genes_anno)
dim(genes_anno)
head(genes_anno)
class(genes_anno)
write.table(genes_anno, file="genes_anno_GSE21034.txt", sep="\t", quote=F, row.names=F)

genes_exp <- exp_GSE21034[same,] # expression for selected genes
genes_exp <- as.data.frame(genes_exp)

##Check if row names matched
table(row.names(genes_exp) == row.names(genes_anno))



# Convert row names into gene names
genes_exp$geneName <- genes_anno$gene_name
#View(as.data.frame(genes_exp[ , 548]))
rownames(genes_exp) <- genes_exp$geneName
genes_exp <- genes_exp[ , -ncol(genes_exp)]
View(genes_exp)
dim(genes_exp)
head(genes_exp)
class(genes_exp)
##


# Transpose Expression matrix
data <- t(genes_exp)
View(data)
class(data)


# Extract clinical information
class(pheno_GSE21034)
View(pheno_GSE21034)
Type <- pheno_GSE21034[ , colnames(pheno_GSE21034) %in% c("sample_id", "pcadb_group")] 
View(Type)

# same_type=intersect(row.names(data),row.names(Type)) #以'df_clinic.genemarkers' 的行名来找因为它在前而且他是最一开始的行名顺序，与‘combind_model...’一致的
# length(same_type)
# df1=data[same_type,] # expression for selected genes
# View(df1)
# dim(df1)
# head(df1)
# class(df1)
# 
# df2 <- Type[same_type,] # pheno of the samples
# View(df2)
# class(df2)
# 
# ##Check if row names matched
# table(row.names(df1) == row.names(df2))

exp_GSE21034_genes <- cbind(data,Type)
dim(exp_GSE21034_genes)
View(exp_GSE21034_genes)


#输出目标基因的表达量
outTab=exp_GSE21034_genes
names(outTab)[ncol(outTab)] <- "Type"
View(outTab)
outTab=cbind(ID=row.names(outTab), outTab)
write.table(outTab, file="geneExp_GSE21034.txt", sep="\t", quote=F, row.names=F)

# This dataset has multiple sample types, but we only want primary and normal
table(outTab$Type)
outTab <- outTab[ outTab$Type %in% c("Primary", "Normal"), ] # filtering for primary and normal only


# Differential analysis


# LOOP for diff analysis

genes <- rownames(genes_exp)

for(gene in genes) {
  
  exp <- outTab [ , c(gene,"Type")] 
  View(exp)
  str(exp)
  group=levels(factor(exp$Type))
  exp$Type=factor(exp$Type, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  #create boxplot
  boxplot=ggboxplot(exp, x="Type", y=gene, color="Type",
                    xlab="",
                    ylab=paste0(gene, " expression"),
                    legend.title="Type",
                    palette = c("blue","red"),
                    add = "jitter")+ 
    stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
  
  #save the plot
  pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
  print(boxplot)
  dev.off()
  
}









# ## "ZWINT"
# gene <- "ZWINT"
# exp <- outTab [ , c("ZWINT","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="ZWINT", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ## "PCBP1"
# 
# gene <- "PCBP1"
# exp <- outTab [ , c("PCBP1","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="PCBP1", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# ## "PABPN1" 
# 
# gene <- "PABPN1"
# exp <- outTab [ , c("PABPN1","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="PABPN1", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ## "PTPRF"
# gene <- "PTPRF"
# exp <- outTab [ , c("PTPRF","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="PTPRF", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# 
# ## "DANCR"
# 
# gene <- "DANCR"
# exp <- outTab [ , c("DANCR","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="DANCR", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# ## "MYC"
# 
# 
# gene <- "MYC"
# exp <- outTab [ , c("MYC","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="MYC", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ## "EEF1D"
# 
# 
# gene <- "EEF1D"
# exp <- outTab [ , c("EEF1D","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="EEF1D", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# ## "FBP1"
# 
# gene <- "FBP1"
# exp <- outTab [ , c("FBP1","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="FBP1", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ##"LAMTOR2"
# 
# 
# 
# gene <- "LAMTOR2"
# exp <- outTab [ , c("LAMTOR2","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="LAMTOR2", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()







##################################################################################################################################################################################
#################################################################### GSE70768 (Cambridge) ########################################################################################
##################################################################################################################################################################################

GSE70768 <- readRDS('Cambridge_eSet.RDS')
GSE70768
View(exprs(GSE70768))
View(pData(GSE70768))
View(fData(GSE70768))
boxplot(GSE70768)



# Expression
exp_GSE70768 <- exprs(GSE70768)
View(exp_GSE70768)
dim(exp_GSE70768)
write.table(exp_GSE70768, file="exp_GSE70768.txt", sep="\t", quote=F, row.names=F)

# Clinical
clinic_GSE70768 <- pData(GSE70768)
dim(clinic_GSE70768)
write.table(clinic_GSE70768, file="target_GSE70768.txt", sep="\t", quote=F, row.names=F)


# Annotation
anno <- readRDS("PCaDB_Gene_Annotation.RDS")
anno_spec <- anno[ anno$gene_name %in% c("ZWINT", 
                                         "APEX1",
                                         "APOF",
                                         "BZW2",
                                         "CACNA1D",
                                         "CAMKK2",
                                         "DANCR",
                                         "DDT",
                                         "DDTL",
                                         "DPM2",
                                         "EEF1D",
                                         "ERBB3",
                                         "FASN",
                                         "FBP1",
                                         "FKBP4",
                                         "GJB1",
                                         "GLT8D2",
                                         "GMDS",
                                         "HOXC6",
                                         "IMPDH2",
                                         "LAMTOR2",
                                         "MRPL24",
                                         "MYC",
                                         "MYL5",
                                         "NOP16",
                                         "OR51F2",
                                         "PGM5P3-AS1",
                                         "PHB2",
                                         "PILRB",
                                         "POLD2",
                                         "PSMG4",
                                         "PTPRF",
                                         "RAB25",
                                         "RPLP0",
                                         "SDF4",
                                         "SLC7A11",
                                         "SMIM22",
                                         "SMUG1",
                                         "SNHG19",
                                         "SNORD49A",
                                         "TMED3",
                                         "TMEM183A",
                                         "TMEM97",
                                         "TPI1",
                                         "TRPM4",
                                         "TRPV6",
                                         "TTLL12",
                                         "ZNHIT1"), ] # for genes of interest only



View(anno_spec)
dim(anno_spec) #48*27



# Processing for single genes
pheno_GSE70768 = read.table("target_GSE70768.txt",header = T,sep="\t")
rownames(pheno_GSE70768) <- pheno_GSE70768$sample_id
View(pheno_GSE70768)  
unique(pheno_GSE70768$pcadb_group)
table(pheno_GSE70768$pcadb_group)
table(colnames(GSE21034) == rownames(pheno_GSE21034)) # Note: the PCaDB annotated the sample id differently, wheras the original expression data has the sample id with GSMxxx from the GEO



# Assigning gene names to the rows of expression matrix
same=intersect(row.names(anno_spec),row.names(exp_GSE70768)) #以'df_clinic.genemarkers' 的行名来找因为它在前而且他是最一开始的行名顺序，与‘combind_model...’一致的

genes_anno=anno_spec[same,] # annotation for selected genes
View(genes_anno)
dim(genes_anno)
head(genes_anno)
class(genes_anno)
write.table(genes_anno, file="genes_anno_GSE70768.txt", sep="\t", quote=F, row.names=F)

genes_exp <- exp_GSE70768[same,] # expression for selected genes
genes_exp <- as.data.frame(genes_exp)

##Check if row names matched
table(row.names(genes_exp) == row.names(genes_anno))



# Convert row names into gene names
genes_exp$geneName <- genes_anno$gene_name
#View(as.data.frame(genes_exp[ , 548]))
rownames(genes_exp) <- genes_exp$geneName
genes_exp <- genes_exp[ , -ncol(genes_exp)]
View(genes_exp)
dim(genes_exp)
head(genes_exp)
class(genes_exp)
##


# Transpose Expression matrix
data <- t(genes_exp)
View(data)
class(data)


# Extract clinical information
class(pheno_GSE70768)
View(pheno_GSE70768)
Type <- pheno_GSE70768[ , colnames(pheno_GSE70768) %in% c("sample_id", "pcadb_group")] 
View(Type)

# same_type=intersect(row.names(data),row.names(Type)) #以'df_clinic.genemarkers' 的行名来找因为它在前而且他是最一开始的行名顺序，与‘combind_model...’一致的
# length(same_type)
# df1=data[same_type,] # expression for selected genes
# View(df1)
# dim(df1)
# head(df1)
# class(df1)
# 
# df2 <- Type[same_type,] # pheno of the samples
# View(df2)
# class(df2)
# 
# ##Check if row names matched
# table(row.names(df1) == row.names(df2))

exp_GSE70768_genes <- cbind(data,Type)
dim(exp_GSE70768_genes)
View(exp_GSE70768_genes)


#输出目标基因的表达量
outTab=exp_GSE70768_genes
names(outTab)[ncol(outTab)] <- "Type"
View(outTab)
outTab=cbind(ID=row.names(outTab), outTab)
write.table(outTab, file="geneExp_GSE70768.txt", sep="\t", quote=F, row.names=F)

# This dataset has multiple sample types, but we only want primary and normal
table(outTab$Type)
outTab <- outTab[ outTab$Type %in% c("Primary", "Normal"), ] # filtering for primary and normal only




# Differential analysis


# LOOP for diff analysis

genes <- rownames(genes_exp)

for(gene in genes) {
  
  exp <- outTab [ , c(gene,"Type")] 
  View(exp)
  str(exp)
  group=levels(factor(exp$Type))
  exp$Type=factor(exp$Type, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  #create boxplot
  boxplot=ggboxplot(exp, x="Type", y=gene, color="Type",
                    xlab="",
                    ylab=paste0(gene, " expression"),
                    legend.title="Type",
                    palette = c("blue","red"),
                    add = "jitter")+ 
    stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
  
  #save the plot
  pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
  print(boxplot)
  dev.off()
  
}

# 
# 
# ## "ZWINT"
# gene <- "ZWINT"
# exp <- outTab [ , c("ZWINT","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="ZWINT", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ## "PCBP1"
# 
# gene <- "PCBP1"
# exp <- outTab [ , c("PCBP1","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="PCBP1", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# ## "PABPN1" 
# 
# gene <- "PABPN1"
# exp <- outTab [ , c("PABPN1","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="PABPN1", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ## "PTPRF"
# gene <- "PTPRF"
# exp <- outTab [ , c("PTPRF","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="PTPRF", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# 
# ## "DANCR"
# 
# gene <- "DANCR"
# exp <- outTab [ , c("DANCR","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="DANCR", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# ## "MYC"
# 
# 
# gene <- "MYC"
# exp <- outTab [ , c("MYC","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="MYC", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ## "EEF1D"
# 
# 
# gene <- "EEF1D"
# exp <- outTab [ , c("EEF1D","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="EEF1D", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# ## "FBP1"
# 
# gene <- "FBP1"
# exp <- outTab [ , c("FBP1","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="FBP1", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ##"LAMTOR2"
# 
# 
# 
# gene <- "LAMTOR2"
# exp <- outTab [ , c("LAMTOR2","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="LAMTOR2", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()








##################################################################################################################################################################################
################################################################### GSE94767 (CancerMap) #########################################################################################
##################################################################################################################################################################################
GSE94767 <- readRDS('CancerMap_eSet.RDS')
GSE94767
View(exprs(GSE94767))
View(pData(GSE94767))
View(fData(GSE94767))
boxplot(GSE94767)



# Expression
exp_GSE94767 <- exprs(GSE94767)
View(exp_GSE94767)
dim(exp_GSE94767)
write.table(exp_GSE94767, file="exp_GSE94767.txt", sep="\t", quote=F, row.names=F)

# Clinical
clinic_GSE94767 <- pData(GSE94767)
dim(clinic_GSE94767)
write.table(clinic_GSE94767, file="target_GSE94767.txt", sep="\t", quote=F, row.names=F)


# Annotation
anno <- readRDS("PCaDB_Gene_Annotation.RDS")
anno_spec <- anno[ anno$gene_name %in% c("ZWINT", 
                                         "APEX1",
                                         "APOF",
                                         "BZW2",
                                         "CACNA1D",
                                         "CAMKK2",
                                         "DANCR",
                                         "DDT",
                                         "DDTL",
                                         "DPM2",
                                         "EEF1D",
                                         "ERBB3",
                                         "FASN",
                                         "FBP1",
                                         "FKBP4",
                                         "GJB1",
                                         "GLT8D2",
                                         "GMDS",
                                         "HOXC6",
                                         "IMPDH2",
                                         "LAMTOR2",
                                         "MRPL24",
                                         "MYC",
                                         "MYL5",
                                         "NOP16",
                                         "OR51F2",
                                         "PGM5P3-AS1",
                                         "PHB2",
                                         "PILRB",
                                         "POLD2",
                                         "PSMG4",
                                         "PTPRF",
                                         "RAB25",
                                         "RPLP0",
                                         "SDF4",
                                         "SLC7A11",
                                         "SMIM22",
                                         "SMUG1",
                                         "SNHG19",
                                         "SNORD49A",
                                         "TMED3",
                                         "TMEM183A",
                                         "TMEM97",
                                         "TPI1",
                                         "TRPM4",
                                         "TRPV6",
                                         "TTLL12",
                                         "ZNHIT1"), ] # for genes of interest only


View(anno_spec)
dim(anno_spec) #48*27



# Processing for single genes
pheno_GSE94767 = read.table("target_GSE94767.txt",header = T,sep="\t")
rownames(pheno_GSE94767) <- pheno_GSE94767$sample_id
pheno_GSE94767$sample_id[is.na(pheno_GSE94767$sample_id)] <- "Unknown" # there is a miising sample so assign it a name
View(pheno_GSE94767)  
unique(pheno_GSE94767$pcadb_group)
table(pheno_GSE94767$pcadb_group)
table(colnames(GSE94767) == rownames(pheno_GSE94767)) # Note: the PCaDB annotated the sample id differently, wheras the original expression data has the sample id with GSMxxx from the GEO



# Assigning gene names to the rows of expression matrix
same=intersect(row.names(anno_spec),row.names(exp_GSE94767)) #以'df_clinic.genemarkers' 的行名来找因为它在前而且他是最一开始的行名顺序，与‘combind_model...’一致的

genes_anno=anno_spec[same,] # annotation for selected genes
View(genes_anno)
dim(genes_anno)
head(genes_anno)
class(genes_anno)

write.table(genes_anno, file="genes_anno_GSE94767.txt", sep="\t", quote=F, row.names=F)





genes_exp <- exp_GSE94767[same,] # expression for selected genes
genes_exp <- as.data.frame(genes_exp)

##Check if row names matched
table(row.names(genes_exp) == row.names(genes_anno)) # Note: in GSE94767, DANCR was not detected



# Convert row names into gene names
genes_exp$geneName <- genes_anno$gene_name
#View(as.data.frame(genes_exp[ , 548]))
rownames(genes_exp) <- genes_exp$geneName
genes_exp <- genes_exp[ , -ncol(genes_exp)]
View(genes_exp)
dim(genes_exp)
head(genes_exp)
class(genes_exp)
##


# Transpose Expression matrix
data <- t(genes_exp)
View(data)
class(data)


# Extract clinical information
class(pheno_GSE94767)
View(pheno_GSE94767)
Type <- pheno_GSE94767[ , colnames(pheno_GSE94767) %in% c("sample_id", "pcadb_group")] 
View(Type)

# same_type=intersect(row.names(data),row.names(Type)) #以'df_clinic.genemarkers' 的行名来找因为它在前而且他是最一开始的行名顺序，与‘combind_model...’一致的
# length(same_type)
# df1=data[same_type,] # expression for selected genes
# View(df1)
# dim(df1)
# head(df1)
# class(df1)
# 
# df2 <- Type[same_type,] # pheno of the samples
# View(df2)
# class(df2)
# 
# ##Check if row names matched
# table(row.names(df1) == row.names(df2))

exp_GSE94767_genes <- cbind(data,Type)
dim(exp_GSE94767_genes)
View(exp_GSE94767_genes)


#输出目标基因的表达量
outTab=exp_GSE94767_genes
names(outTab)[ncol(outTab)] <- "Type"
View(outTab)
outTab=cbind(ID=row.names(outTab), outTab)
write.table(outTab, file="geneExp_GSE94767.txt", sep="\t", quote=F, row.names=F)

# This dataset has multiple sample types, but we only want primary and normal
table(outTab$Type)
outTab <- outTab[ outTab$Type %in% c("Tumor", "Normal"), ] # filtering for tumour and normal only
table(outTab$Type)


# Differential analysis


# LOOP for diff analysis

genes <- rownames(genes_exp)

for(gene in genes) {
  
  exp <- outTab [ , c(gene,"Type")] 
  View(exp)
  str(exp)
  group=levels(factor(exp$Type))
  exp$Type=factor(exp$Type, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  #create boxplot
  boxplot=ggboxplot(exp, x="Type", y=gene, color="Type",
                    xlab="",
                    ylab=paste0(gene, " expression"),
                    legend.title="Type",
                    palette = c("blue","red"),
                    add = "jitter")+ 
    stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
  
  #save the plot
  pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
  print(boxplot)
  dev.off()
  
}





# ## "ZWINT"
# gene <- "ZWINT"
# exp <- outTab [ , c("ZWINT","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="ZWINT", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ## "PCBP1"
# 
# gene <- "PCBP1"
# exp <- outTab [ , c("PCBP1","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="PCBP1", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# ## "PABPN1" 
# 
# gene <- "PABPN1"
# exp <- outTab [ , c("PABPN1","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="PABPN1", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ## "PTPRF"
# gene <- "PTPRF"
# exp <- outTab [ , c("PTPRF","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="PTPRF", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# 
# # ## "DANCR"
# # 
# # gene <- "DANCR"
# # exp <- outTab [ , c("DANCR","Type")] 
# # View(exp)
# # str(exp)
# # group=levels(factor(exp$Type))
# # exp$Type=factor(exp$Type, levels=group)
# # comp=combn(group,2)
# # my_comparisons=list()
# # for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# # 
# # #绘制boxplot
# # boxplot=ggboxplot(exp, x="Type", y="DANCR", color="Type",
# #                   xlab="",
# #                   ylab=paste0(gene, " expression"),
# #                   legend.title="Type",
# #                   palette = c("blue","red"),
# #                   add = "jitter")+ 
# #   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# # 
# # #输出图片
# # pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# # print(boxplot)
# # dev.off()
# 
# 
# ## "MYC"
# 
# 
# gene <- "MYC"
# exp <- outTab [ , c("MYC","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="MYC", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ## "EEF1D"
# 
# 
# gene <- "EEF1D"
# exp <- outTab [ , c("EEF1D","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="EEF1D", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# ## "FBP1"
# 
# gene <- "FBP1"
# exp <- outTab [ , c("FBP1","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="FBP1", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ##"LAMTOR2"
# 
# 
# 
# gene <- "LAMTOR2"
# exp <- outTab [ , c("LAMTOR2","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="LAMTOR2", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()






##################################################################################################################################################################################
################################################################## E-MTAB-6128 (CIT) #############################################################################################
##################################################################################################################################################################################


CIT <- readRDS('CIT_eSet.RDS')
CIT
View(exprs(CIT))
View(pData(CIT))
View(fData(CIT))
boxplot(CIT)



# Expression
exp_CIT <- exprs(CIT)
View(exp_CIT)
dim(exp_CIT)
write.table(exp_CIT, file="exp_CIT.txt", sep="\t", quote=F, row.names=F)

# Clinical
clinic_CIT <- pData(CIT)
dim(clinic_CIT)
write.table(clinic_CIT, file="target_CIT.txt", sep="\t", quote=F, row.names=F)


# Annotation
anno <- readRDS("PCaDB_Gene_Annotation.RDS")
anno_spec <- anno[ anno$gene_name %in% c("ZWINT", 
                                         "APEX1",
                                         "APOF",
                                         "BZW2",
                                         "CACNA1D",
                                         "CAMKK2",
                                         "DANCR",
                                         "DDT",
                                         "DDTL",
                                         "DPM2",
                                         "EEF1D",
                                         "ERBB3",
                                         "FASN",
                                         "FBP1",
                                         "FKBP4",
                                         "GJB1",
                                         "GLT8D2",
                                         "GMDS",
                                         "HOXC6",
                                         "IMPDH2",
                                         "LAMTOR2",
                                         "MRPL24",
                                         "MYC",
                                         "MYL5",
                                         "NOP16",
                                         "OR51F2",
                                         "PGM5P3-AS1",
                                         "PHB2",
                                         "PILRB",
                                         "POLD2",
                                         "PSMG4",
                                         "PTPRF",
                                         "RAB25",
                                         "RPLP0",
                                         "SDF4",
                                         "SLC7A11",
                                         "SMIM22",
                                         "SMUG1",
                                         "SNHG19",
                                         "SNORD49A",
                                         "TMED3",
                                         "TMEM183A",
                                         "TMEM97",
                                         "TPI1",
                                         "TRPM4",
                                         "TRPV6",
                                         "TTLL12",
                                         "ZNHIT1"), ] # for genes of interest only




View(anno_spec)
dim(anno_spec) #48*27



# Processing for single genes
pheno_CIT = read.table("target_CIT.txt",header = T,sep="\t")
rownames(pheno_CIT) <- pheno_CIT$sample_id
View(pheno_CIT)  
unique(pheno_CIT$pcadb_group)
table(pheno_CIT$pcadb_group)
table(colnames(CIT) == rownames(pheno_CIT)) # Note: the PCaDB annotated the sample id differently, wheras the original expression data has the sample id with GSMxxx from the GEO



# Assigning gene names to the rows of expression matrix
same=intersect(row.names(anno_spec),row.names(exp_CIT)) #以'df_clinic.genemarkers' 的行名来找因为它在前而且他是最一开始的行名顺序，与‘combind_model...’一致的

genes_anno=anno_spec[same,] # annotation for selected genes
View(genes_anno)
dim(genes_anno)
head(genes_anno)
class(genes_anno)
write.table(genes_anno, file="genes_anno_CIT.txt", sep="\t", quote=F, row.names=F)
genes_exp <- exp_CIT[same,] # expression for selected genes
genes_exp <- as.data.frame(genes_exp)

##Check if row names matched
table(row.names(genes_exp) == row.names(genes_anno))



# Convert row names into gene names
genes_exp$geneName <- genes_anno$gene_name
#View(as.data.frame(genes_exp[ , 548]))
rownames(genes_exp) <- genes_exp$geneName
genes_exp <- genes_exp[ , -ncol(genes_exp)]
View(genes_exp)
dim(genes_exp)
head(genes_exp)
class(genes_exp)
##


# Transpose Expression matrix
data <- t(genes_exp)
View(data)
class(data)


# Extract clinical information
class(pheno_CIT)
View(pheno_CIT)
Type <- pheno_CIT[ , colnames(pheno_CIT) %in% c("sample_id", "pcadb_group")] 
View(Type)

# same_type=intersect(row.names(data),row.names(Type)) #以'df_clinic.genemarkers' 的行名来找因为它在前而且他是最一开始的行名顺序，与‘combind_model...’一致的
# length(same_type)
# df1=data[same_type,] # expression for selected genes
# View(df1)
# dim(df1)
# head(df1)
# class(df1)
# 
# df2 <- Type[same_type,] # pheno of the samples
# View(df2)
# class(df2)
# 
# ##Check if row names matched
# table(row.names(df1) == row.names(df2))

exp_CIT_genes <- cbind(data,Type)
dim(exp_CIT_genes)
View(exp_CIT_genes)


#输出目标基因的表达量
outTab=exp_CIT_genes
names(outTab)[ncol(outTab)] <- "Type"
View(outTab)
outTab=cbind(ID=row.names(outTab), outTab)
write.table(outTab, file="geneExp_CIT.txt", sep="\t", quote=F, row.names=F)

# This dataset has multiple sample types, but we only want primary and normal
table(outTab$Type)
outTab <- outTab[ outTab$Type %in% c("Primary", "Normal"), ] # filtering for primary and normal only





# Differential analysis



# LOOP for diff analysis

genes <- rownames(genes_exp)

for(gene in genes) {
  
  exp <- outTab [ , c(gene,"Type")] 
  View(exp)
  str(exp)
  group=levels(factor(exp$Type))
  exp$Type=factor(exp$Type, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  
  #create boxplot
  boxplot=ggboxplot(exp, x="Type", y=gene, color="Type",
                    xlab="",
                    ylab=paste0(gene, " expression"),
                    legend.title="Type",
                    palette = c("blue","red"),
                    add = "jitter")+ 
    stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
  
  #save the plot
  pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
  print(boxplot)
  dev.off()
  
}


# ## "ZWINT"
# gene <- "ZWINT"
# exp <- outTab [ , c("ZWINT","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="ZWINT", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ## "PCBP1"
# 
# gene <- "PCBP1"
# exp <- outTab [ , c("PCBP1","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="PCBP1", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# ## "PABPN1" 
# 
# gene <- "PABPN1"
# exp <- outTab [ , c("PABPN1","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="PABPN1", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ## "PTPRF"
# gene <- "PTPRF"
# exp <- outTab [ , c("PTPRF","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="PTPRF", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# 
# ## "DANCR"
# 
# gene <- "DANCR"
# exp <- outTab [ , c("DANCR","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="DANCR", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# ## "MYC"
# 
# 
# gene <- "MYC"
# exp <- outTab [ , c("MYC","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="MYC", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ## "EEF1D"
# 
# 
# gene <- "EEF1D"
# exp <- outTab [ , c("EEF1D","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="EEF1D", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# ## "FBP1"
# 
# gene <- "FBP1"
# exp <- outTab [ , c("FBP1","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="FBP1", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()
# 
# 
# 
# ##"LAMTOR2"
# 
# 
# 
# gene <- "LAMTOR2"
# exp <- outTab [ , c("LAMTOR2","Type")] 
# View(exp)
# str(exp)
# group=levels(factor(exp$Type))
# exp$Type=factor(exp$Type, levels=group)
# comp=combn(group,2)
# my_comparisons=list()
# for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
# 
# #绘制boxplot
# boxplot=ggboxplot(exp, x="Type", y="LAMTOR2", color="Type",
#                   xlab="",
#                   ylab=paste0(gene, " expression"),
#                   legend.title="Type",
#                   palette = c("blue","red"),
#                   add = "jitter")+ 
#   stat_compare_means(comparisons=my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif", method = "wilcox.test")
# 
# #输出图片
# pdf(file=paste0(gene,".diff.pdf"), width=5, height=4.5)
# print(boxplot)
# dev.off()





#-----------------------------------------------------------------##
#                          GeneMeta                               ##
#-----------------------------------------------------------------##   
library(GeneMeta)
library(RColorBrewer)
library(Biobase)
library(readxl)
library(stringr)
library(impute)


##################################################################################################################################################################################
################################################################## For all  genes ###############################################################################################
##################################################################################################################################################################################

## TCGA
data_TCGA = read.table("/Users/zhuofanmou/Documents/Phd/papers/paper3/data/GeneMeta_prostate_tissue/meta_analysis/differential_v2/TCGA_res/geneExp_TCGA.txt",header = T,sep="\t",row.names = 1)
View((data_TCGA))

table(data_TCGA$Type)
data_TCGA <- data_TCGA[ data_TCGA$Type %in% c("Primary","Normal") , ]
View((data_TCGA))

temp <- data_TCGA

data_TCGA <- data_TCGA[ , -c(ncol(data_TCGA), ncol(data_TCGA) - 1)]
data_TCGA <-t(data_TCGA)
data_TCGA <- as.data.frame(data_TCGA)
class(data_TCGA)
View((data_TCGA))

# read in phenotype
pheno0_TCGA = read.table("/Users/zhuofanmou/Documents/Phd/papers/paper3/data/GeneMeta_prostate_tissue/meta_analysis/differential_v2/TCGA_res/geneExp_TCGA.txt",header = T,sep="\t")
View(head(pheno0_TCGA))
pheno0_TCGA <- pheno0_TCGA[ , c(1, ncol(pheno0_TCGA))]
colnames(pheno0_TCGA) <- c("Sample","Group")
pheno0_TCGA[,2] = str_remove_all(pheno0_TCGA[,2]," +$")

pheno_TCGA = pheno0_TCGA[match(colnames(data_TCGA),pheno0_TCGA[,1]),]
rownames(pheno_TCGA) <-pheno_TCGA$Sample
View((pheno_TCGA))

same=intersect(row.names(temp),row.names(pheno_TCGA)) #以'df_clinic.genemarkers' 的行名来找因为它在前而且他是最一开始的行名顺序，与‘combind_model...’一致的
length(same)
pheno_TCGA <- pheno_TCGA[same,] # pheno of the samples
View(pheno_TCGA)
class(pheno_TCGA)

##Check if row names matched
table(row.names(pheno_TCGA) == row.names(temp))



# GSE21034


data_GSE21034 = read.table("/Users/zhuofanmou/Documents/Phd/papers/paper3/data/GeneMeta_prostate_tissue/meta_analysis/differential_v2/GSE21034_res/geneExp_GSE21034.txt",header = T,sep="\t",row.names = 1)
View((data_GSE21034))

table(data_GSE21034$Type)
data_GSE21034 <- data_GSE21034[ data_GSE21034$Type %in% c("Primary","Normal") , ]
View((data_GSE21034))

temp <- data_GSE21034

data_GSE21034 <- data_GSE21034[ , -c(ncol(data_GSE21034), ncol(data_GSE21034) - 1)]
data_GSE21034 <-t(data_GSE21034)
data_GSE21034 <- as.data.frame(data_GSE21034)
class(data_GSE21034)
View((data_GSE21034))

# read in phenotype
pheno0_GSE21034 = read.table("/Users/zhuofanmou/Documents/Phd/papers/paper3/data/GeneMeta_prostate_tissue/meta_analysis/differential_v2/GSE21034_res/geneExp_GSE21034.txt",header = T,sep="\t")
View(head(pheno0_GSE21034))
pheno0_GSE21034 <- pheno0_GSE21034[ , c(1, ncol(pheno0_GSE21034))]
colnames(pheno0_GSE21034) <- c("Sample","Group")
pheno0_GSE21034[,2] = str_remove_all(pheno0_GSE21034[,2]," +$")

pheno_GSE21034 = pheno0_GSE21034[match(colnames(data_GSE21034),pheno0_GSE21034[,1]),]
rownames(pheno_GSE21034) <-pheno_GSE21034$Sample
View((pheno_GSE21034))


same=intersect(row.names(temp),row.names(pheno_GSE21034)) #以'df_clinic.genemarkers' 的行名来找因为它在前而且他是最一开始的行名顺序，与‘combind_model...’一致的

length(same)
pheno_GSE21034 <- pheno_GSE21034[same,] # pheno of the samples
View(pheno_GSE21034)
class(pheno_GSE21034)

##Check if row names matched
table(row.names(pheno_GSE21034) == row.names(temp))



# GSE70768

data_GSE70768 = read.table("/Users/zhuofanmou/Documents/Phd/papers/paper3/data/GeneMeta_prostate_tissue/meta_analysis/differential_v2/GSE70768_res/geneExp_GSE70768.txt",header = T,sep="\t",row.names = 1)
View((data_GSE70768))

table(data_GSE70768$Type)
data_GSE70768 <- data_GSE70768[ data_GSE70768$Type %in% c("Primary","Normal") , ]
View((data_GSE70768))

temp <- data_GSE70768

data_GSE70768 <- data_GSE70768[ , -c(ncol(data_GSE70768), ncol(data_GSE70768) - 1)]
data_GSE70768 <-t(data_GSE70768)
data_GSE70768 <- as.data.frame(data_GSE70768)
class(data_GSE70768)
View((data_GSE70768))

# read in phenotype
pheno0_GSE70768 = read.table("/Users/zhuofanmou/Documents/Phd/papers/paper3/data/GeneMeta_prostate_tissue/meta_analysis/differential_v2/GSE70768_res/geneExp_GSE70768.txt",header = T,sep="\t")
View(head(pheno0_GSE70768))
pheno0_GSE70768 <- pheno0_GSE70768[ , c(1, ncol(pheno0_GSE70768))]
colnames(pheno0_GSE70768) <- c("Sample","Group")
pheno0_GSE70768[,2] = str_remove_all(pheno0_GSE70768[,2]," +$")

pheno_GSE70768 = pheno0_GSE70768[match(colnames(data_GSE70768),pheno0_GSE70768[,1]),]
rownames(pheno_GSE70768) <-pheno_GSE70768$Sample
View((pheno_GSE70768))


same=intersect(row.names(temp),row.names(pheno_GSE70768)) #以'df_clinic.genemarkers' 的行名来找因为它在前而且他是最一开始的行名顺序，与‘combind_model...’一致的

length(same)
pheno_GSE70768 <- pheno_GSE70768[same,] # pheno of the samples
View(pheno_GSE70768)
class(pheno_GSE70768)

##Check if row names matched
table(row.names(pheno_GSE70768) == row.names(temp))







# GSE94767 

data_GSE94767 = read.table("/Users/zhuofanmou/Documents/Phd/papers/paper3/data/GeneMeta_prostate_tissue/meta_analysis/differential_v2/GSE94767_res/geneExp_GSE94767.txt",header = T,sep="\t",row.names = 1)
View((data_GSE94767))

table(data_GSE94767$Type)
data_GSE94767 <- data_GSE94767[ data_GSE94767$Type %in% c("Tumor","Normal") , ]
View((data_GSE94767))

temp <- data_GSE94767

data_GSE94767 <- data_GSE94767[ , -c(ncol(data_GSE94767), ncol(data_GSE94767) - 1)]
data_GSE94767 <-t(data_GSE94767)
data_GSE94767 <- as.data.frame(data_GSE94767)
class(data_GSE94767)
View((data_GSE94767))

# read in phenotype
pheno0_GSE94767 = read.table("/Users/zhuofanmou/Documents/Phd/papers/paper3/data/GeneMeta_prostate_tissue/meta_analysis/differential_v2/GSE94767_res/geneExp_GSE94767.txt",header = T,sep="\t")
View(head(pheno0_GSE94767))
pheno0_GSE94767 <- pheno0_GSE94767[ , c(1, ncol(pheno0_GSE94767))]
colnames(pheno0_GSE94767) <- c("Sample","Group")
pheno0_GSE94767[,2] = str_remove_all(pheno0_GSE94767[,2]," +$")

pheno_GSE94767 = pheno0_GSE94767[match(colnames(data_GSE94767),pheno0_GSE94767[,1]),]
rownames(pheno_GSE94767) <-pheno_GSE94767$Sample
View((pheno_GSE94767))


same=intersect(row.names(temp),row.names(pheno_GSE94767)) #以'df_clinic.genemarkers' 的行名来找因为它在前而且他是最一开始的行名顺序，与‘combind_model...’一致的

length(same)
pheno_GSE94767 <- pheno_GSE94767[same,] # pheno of the samples
View(pheno_GSE94767)
class(pheno_GSE94767)

##Check if row names matched
table(row.names(pheno_GSE94767) == row.names(temp))



# CIT

data_CIT = read.table("/Users/zhuofanmou/Documents/Phd/papers/paper3/data/GeneMeta_prostate_tissue/meta_analysis/differential_v2/CIT_res/geneExp_CIT.txt",header = T,sep="\t",row.names = 1)
View((data_CIT))

table(data_CIT$Type)
data_CIT <- data_CIT[ data_CIT$Type %in% c("Primary","Normal") , ]
View((data_CIT))

temp <- data_CIT

data_CIT <- data_CIT[ , -c(ncol(data_CIT), ncol(data_CIT) - 1)]
data_CIT <-t(data_CIT)
data_CIT <- as.data.frame(data_CIT)
class(data_CIT)
View((data_CIT))

# read in phenotype
pheno0_CIT = read.table("/Users/zhuofanmou/Documents/Phd/papers/paper3/data/GeneMeta_prostate_tissue/meta_analysis/differential_v2/CIT_res/geneExp_CIT.txt",header = T,sep="\t")
View(head(pheno0_CIT))
pheno0_CIT <- pheno0_CIT[ , c(1, ncol(pheno0_CIT))]
colnames(pheno0_CIT) <- c("Sample","Group")
pheno0_CIT[,2] = str_remove_all(pheno0_CIT[,2]," +$")

pheno_CIT = pheno0_CIT[match(colnames(data_CIT),pheno0_CIT[,1]),]
rownames(pheno_CIT) <-pheno_CIT$Sample
View((pheno_CIT))


same=intersect(row.names(temp),row.names(pheno_CIT)) #以'df_clinic.genemarkers' 的行名来找因为它在前而且他是最一开始的行名顺序，与‘combind_model...’一致的

length(same)
pheno_CIT <- pheno_CIT[same,] # pheno of the samples
View(pheno_CIT)
class(pheno_CIT)

##Check if row names matched
table(row.names(pheno_CIT) == row.names(temp))



# We need to create ExpressionSet

# TCGA
eset_TCGA <- ExpressionSet(assayData = assayDataNew(exprs = as.matrix(data_TCGA)))
pData(eset_TCGA) <-pheno_TCGA
View(pData(eset_TCGA))
eset_TCGA
eset_TCGA <- eset_TCGA[!(rownames(eset_TCGA)=="DANCR") , ]
eset_TCGA <- eset_TCGA[!(rownames(eset_TCGA)=="HOXC6") , ]
eset_TCGA <- eset_TCGA[!(rownames(eset_TCGA)=="OR51F2") , ]
eset_TCGA <- eset_TCGA[!(rownames(eset_TCGA)=="PILRB") , ]
eset_TCGA <- eset_TCGA[!(rownames(eset_TCGA)=="DDT") , ]
eset_TCGA <- eset_TCGA[!(rownames(eset_TCGA)=="SNORD49A") , ]
eset_TCGA <- eset_TCGA[!(rownames(eset_TCGA)=="PGM5P3.AS1") , ]
eset_TCGA <- eset_TCGA[!(rownames(eset_TCGA)=="SNHG19") , ]
eset_TCGA
condition_TCGA <- as.factor(pheno_TCGA$Group)

#levels(condition_TCGA) <- c('Primary', 'Normal')
levels(condition_TCGA) <- c(1,0) #c(0,1) #c(1,0) #c(0,1) # swap so that I think it show the expression changes from tumour to normal; +Mu = up-regulated in tumour, -Mu = down-regulated in tumour
condition_TCGA<- as.numeric(as.character(condition_TCGA))
eset_TCGA

# GSE21034

eset_GSE21034 <- ExpressionSet(assayData = assayDataNew(exprs = as.matrix(data_GSE21034)))
pData(eset_GSE21034) <-pheno_GSE21034
View(pData(eset_GSE21034))
eset_GSE21034
eset_GSE21034 <- eset_GSE21034[!(rownames(eset_GSE21034)=="DANCR") , ]
eset_GSE21034 <- eset_GSE21034[!(rownames(eset_GSE21034)=="HOXC6") , ]
eset_GSE21034 <- eset_GSE21034[!(rownames(eset_GSE21034)=="OR51F2") , ]
eset_GSE21034 <- eset_GSE21034[!(rownames(eset_GSE21034)=="PILRB") , ]
eset_GSE21034 <- eset_GSE21034[!(rownames(eset_GSE21034)=="SNHG19") , ]
eset_GSE21034
condition_GSE21034 <- as.factor(pheno_GSE21034$Group)

#levels(condition_GSE21034) <- c('Primary', 'Normal')
levels(condition_GSE21034) <- c(1,0) #c(0,1)#c(1,0) #c(0,1)
condition_GSE21034<- as.numeric(as.character(condition_GSE21034))
eset_GSE21034

# GSE70768

eset_GSE70768 <- ExpressionSet(assayData = assayDataNew(exprs = as.matrix(data_GSE70768)))
pData(eset_GSE70768) <-pheno_GSE70768
View(pData(eset_GSE70768))
eset_GSE70768
eset_GSE70768 <- eset_GSE70768[!(rownames(eset_GSE70768)=="DANCR") , ]
eset_GSE70768 <- eset_GSE70768[!(rownames(eset_GSE70768)=="HOXC6") , ]
eset_GSE70768 <- eset_GSE70768[!(rownames(eset_GSE70768)=="OR51F2") , ]
eset_GSE70768 <- eset_GSE70768[!(rownames(eset_GSE70768)=="PILRB") , ]
eset_GSE70768 <- eset_GSE70768[!(rownames(eset_GSE70768)=="DDT") , ]
eset_GSE70768 <- eset_GSE70768[!(rownames(eset_GSE70768)=="SNORD49A") , ]
eset_GSE70768
condition_GSE70768 <- as.factor(pheno_GSE70768$Group)

#levels(condition_GSE70768) <- c('Primary', 'Normal')
levels(condition_GSE70768) <- c(1,0) #c(0,1)#c(1,0) #c(0,1)
condition_GSE70768<- as.numeric(as.character(condition_GSE70768))
eset_GSE70768




# GSE94767

eset_GSE94767 <- ExpressionSet(assayData = assayDataNew(exprs = as.matrix(data_GSE94767)))
pData(eset_GSE94767) <-pheno_GSE94767
View(pData(eset_GSE94767))
eset_GSE94767
#eset_GSE94767 <- eset_GSE94767[-(rownames(eset_GSE94767)=="DANCR") , ] # dont need as it does not have DANCR
eset_GSE94767 <- eset_GSE94767[!(rownames(eset_GSE94767)=="SNHG19") , ] 
condition_GSE94767 <- as.factor(pheno_GSE94767$Group)
#levels(condition_GSE94767) <- c('Primary', 'Normal')
levels(condition_GSE94767) <- c(1,0) #c(0,1)#c(1,0) #c(0,1)
condition_GSE94767<- as.numeric(as.character(condition_GSE94767))
eset_GSE94767



# CIT

eset_CIT <- ExpressionSet(assayData = assayDataNew(exprs = as.matrix(data_CIT)))
pData(eset_CIT) <-pheno_CIT
View(pData(eset_CIT))
eset_CIT
eset_CIT <- eset_CIT[!(rownames(eset_CIT)=="DANCR") , ]
eset_CIT <- eset_CIT[!(rownames(eset_CIT)=="HOXC6") , ]
eset_CIT <- eset_CIT[!(rownames(eset_CIT)=="OR51F2") , ]
eset_CIT <- eset_CIT[!(rownames(eset_CIT)=="PILRB") , ]
eset_CIT <- eset_CIT[!(rownames(eset_CIT)=="PGM5P3.AS1") , ]
eset_CIT <- eset_CIT[!(rownames(eset_CIT)=="SNHG19") , ]

condition_CIT <- as.factor(pheno_CIT$Group)
#levels(condition_CIT) <- c('Primary', 'Normal')
levels(condition_CIT) <- c(1,0) #c(0,1) #c(1,0) #c(0,1)
condition_CIT<- as.numeric(as.character(condition_CIT))
eset_CIT



table(featureNames(eset_TCGA) %in% featureNames(eset_CIT))
table(featureNames(eset_TCGA) %in% featureNames(eset_GSE94767))

table(featureNames(eset_TCGA) == featureNames(eset_CIT))
table(featureNames(eset_TCGA) == featureNames(eset_GSE94767))



##################################################################################################################################################################################
################################################################## For all five studies ##########################################################################################
##################################################################################################################################################################################

###################################################
### code chunk number 13: claculationAllinOne
###################################################
esets     <- list(eset_TCGA,
                  eset_GSE21034,
                  eset_GSE70768,
                  eset_GSE94767,
                  eset_CIT)
# data.ER   <-pData(Nevins)[,"ER.status"]
# levels(data.ER) <- c(0,1)
# data.ER<- as.numeric(as.character(data.ER))
classes   <- list(condition_TCGA,
                  condition_GSE21034,
                  condition_GSE70768,
                  condition_GSE94767,
                  condition_CIT)
#theScores <- zScores(esets,classes,useREM=FALSE,CombineExp=1:2)
theScores <- zScores(esets,classes,useREM=TRUE,CombineExp=1:5) # CombineExp can be set to combine wanted datasets
View(theScores)



###################################################
### code chunk number 14: show results
###################################################
View(theScores)
theScores[1:2,]

# output data
#write.table(theScores,"Table1_Meta_synovial.txt",row.names = T,col.names = NA,quote = F,sep="\t")
write.csv(theScores, 'Table1_Meta40.csv')# alternativelt save as csv, probably better
save(theScores, file = "Table1_Meta40.rda") # same the unfiltered normalised expression matrix



# FDR
set.seed(123)
ScoresFDR <- zScoreFDR(esets, classes, useREM=TRUE, nperm=50,CombineExp=1:5)
View(ScoresFDR$two.sided)
#View(ScoresFDR$pos)

ScoresFDR_two.sided <- ScoresFDR$two.sided
View(ScoresFDR_two.sided)
write.csv(ScoresFDR_two.sided, 'Table2_Meta40_FDR_default.csv')# alternativelt save as csv, probably better
save(ScoresFDR_two.sided, file = "Table2_Meta40_FDR_default.rda") # same the unfiltered normalised expression matrix



# BH-FDR

## https://www.youtube.com/watch?v=DRRT4k7A2dM

b<- as.data.frame(ScoresFDR_two.sided)
class(b)
b$rawPbyZ <- 2*pnorm(-abs(b$zSco))
View(b)
b$FDR_BH <- p.adjust(b$rawPbyZ, "BH") # I dont think it is necessary to do fdr if we assuming to looking gene one by one
View(b)
write.csv(b, 'Table3_Meta40_FDR_BH.csv')# alternativelt save as csv, probably better



###################################################
### Pooled ES volcano
###################################################

GeneMeta_res<- as.data.frame(b)
GeneMeta_res$pooled_ES <- GeneMeta_res$MUvals
View(GeneMeta_res)

deg.data <- as.data.frame(GeneMeta_res[ , c(21,20)])
View(deg.data)
class(deg.data)
deg.data$logP <- -log10(deg.data$FDR_BH)

deg.data$Group = "Not significant"
deg.data$Group[which( (deg.data$FDR_BH < 0.05))] = "Significant"
#deg.data$Group[which( (deg.data$adj.P.Val < p.adjust.cutoff) & (deg.data$logFC < -log(fold.change.cutoff,2)) )] = "Down-regulated"
deg.data$Label = ""
deg.data$LabelSF = ""
deg.data <- deg.data[order(deg.data$FDR_BH), ]
up.genes <- head(rownames(deg.data)[which(deg.data$Group == "Significant")], 10)
#down.genes <- head(deg.data$Symbol[which(deg.data$Group == "Down-regulated")], 10)
deg.top10.genes <- (as.character(up.genes))
deg.data$Label[match(deg.top10.genes, rownames(deg.data))] <- deg.top10.genes
deg.data$LabelSF[match(deg.top10.genes, rownames(deg.data))] <- deg.top10.genes # NEED to run ' Extracting splicing factor info' below beforehand for 'degSF'
deg.data$Label_all <- rownames(deg.data)
View(deg.data)

library(ggplot2)
library(ggpubr)
p1 <- ggscatter(deg.data, x = "pooled_ES", y = "logP",
                color = "Group",
                shape = "Group",
                palette = c("#616161", "#000000"),#c("#BBBBBB", "#CC0000"),
                size = 2.5,
                label = deg.data$Label_all,
                #label = deg.data$LabelSF,
                font.label = 10,
                repel = T,
                xlab = "Pooled effect size (Hedges'g)",
                xlim = c(-0.5,1.4),
                ylab = "-log10(BH adjusted p-value)",
                ylim = c(0,55)) +
  theme_bw() +
  geom_hline(yintercept = -log(0.05,10), linetype="dashed") +
  ggtitle("Volcano plot of gene markers")
#geom_vline(xintercept = c(log(fold.change.cutoff,2),-log(fold.change.cutoff,2)), linetype="dashed")+
#ggtitle(paste("Volcano plot","(",colnames(fit2$contrasts),")"))
p1
library(cowplot)
## output plot
save_plot("volcano.pdf",p1,base_width = 12, base_height = 14)
save_plot("volcano.tiff",p1,base_width = 12, base_height = 14,
          compression = "lzw",dpi = 500)

#-----------------------------------------------------------------##
#                          Forest                                 ##
#-----------------------------------------------------------------##   
library(meta)


#https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/pooling-es.html
#https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/effects.html

lname <- load('Table1_Meta40.rda')
lname

DF <- as.data.frame(theScores)
class(DF)
View(DF)


sig_genes <- deg.data [deg.data$FDR_BH < 0.05, ]
View(sig_genes)


# LOOP for meta forest

meta.genes <- rownames(sig_genes)

for(gene in meta.genes) {
  
GEO.study <- c('TCGA','GSE21034','GSE70768', 'GSE94767', 'CIT')

TCGA_DF <- DF[rownames(DF) == gene, ]
TCGA_ES <- TCGA_DF$Effect_Ex_1
TCGA_se <- sqrt(TCGA_DF$EffectVar_Ex_1)
#GSE57218_se <- sqrt(((5+5)/(5*5))+((GSE57218_ES)^2/(2*(5+5))))

GSE21034_DF <- DF[rownames(DF) == gene, ]
GSE21034_ES <- GSE21034_DF$Effect_Ex_2
GSE21034_se <- sqrt(GSE21034_DF$EffectVar_Ex_2)
#GSE117999_se <- sqrt(((10+10)/(10*10))+((GSE117999_ES)^2/(2*(10+10))))


GSE70768_DF <- DF[rownames(DF) == gene, ]
GSE70768_ES <- GSE70768_DF$Effect_Ex_3
GSE70768_se <- sqrt(GSE70768_DF$EffectVar_Ex_3)
#GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))

GSE94767_DF <- DF[rownames(DF) == gene, ]
GSE94767_ES <- GSE94767_DF$Effect_Ex_4
GSE94767_se <- sqrt(GSE94767_DF$EffectVar_Ex_4)
#GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))

CIT_DF <- DF[rownames(DF) == gene, ]
CIT_ES <- CIT_DF$Effect_Ex_5
CIT_se <- sqrt(CIT_DF$EffectVar_Ex_5)
#GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))


## Variable set-up

TE <- c(TCGA_ES, GSE21034_ES , GSE70768_ES, GSE94767_ES, CIT_ES )
seTE <- c(TCGA_se,GSE21034_se, GSE70768_se, GSE94767_se, CIT_se)

# Join the variables to create a data frame
df2_gene <- data.frame(GEO.study,TE,seTE )
df2_gene
rownames(df2_gene) <- df2_gene$GEO.study

data.frame(df2_gene, stringsAsFactors = TRUE)

m.gen <- metagen(TE = TE,
                 seTE = seTE,
                 studlab = GEO.study,
                 data = df2_gene,
                 sm = "SMD",
                 #method.ci = 't',
                 method.smd = "Hedges",
                 fixed = FALSE,
                 random = TRUE,
                 method.tau = "DL",#!!!!!!!!!!! here we use DJ (DerSimonian and Laird estimator) because this is what GeneMeta uses!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 #method.tau = "REML",
                 hakn = FALSE, #if use false then metafor/GeneMeta/meta all matched!!!!!!!!! AND it uses Wald-type test
                 #hakn = TRUE,
                 title = gene)
m.gen
summary(m.gen)
m.gen$pval.random

pdf(file = paste0(gene,".meta_forest.pdf"), width = 14, height = 8)
forest(m.gen, studlab = TRUE)
dev.off()


}





##################################################################################################################################################################################
####################################################### For four studies (TCGA removed) ##########################################################################################
##################################################################################################################################################################################

###################################################
### code chunk number 13: claculationAllinOne
###################################################
esets     <- list(
                  eset_GSE21034,
                  eset_GSE70768,
                  eset_GSE94767,
                  eset_CIT)
# data.ER   <-pData(Nevins)[,"ER.status"]
# levels(data.ER) <- c(0,1)
# data.ER<- as.numeric(as.character(data.ER))
classes   <- list(
                  condition_GSE21034,
                  condition_GSE70768,
                  condition_GSE94767,
                  condition_CIT)
#theScores <- zScores(esets,classes,useREM=FALSE,CombineExp=1:2)
theScores <- zScores(esets,classes,useREM=TRUE,CombineExp=1:4) # CombineExp can be set to combine wanted datasets
View(theScores)



###################################################
### code chunk number 14: show results
###################################################
View(theScores)
theScores[1:2,]

# output data
#write.table(theScores,"Table1_Meta_synovial.txt",row.names = T,col.names = NA,quote = F,sep="\t")
write.csv(theScores, 'Table1_Meta40_TCGArem.csv')# alternativelt save as csv, probably better
save(theScores, file = "Table1_Meta40_TCGArem.rda") # same the unfiltered normalised expression matrix



# FDR
set.seed(123)
ScoresFDR <- zScoreFDR(esets, classes, useREM=TRUE, nperm=50,CombineExp=1:4)
View(ScoresFDR$two.sided)
#View(ScoresFDR$pos)

ScoresFDR_two.sided <- ScoresFDR$two.sided
View(ScoresFDR_two.sided)
write.csv(ScoresFDR_two.sided, 'Table2_Meta40_FDR_default_TCGArem.csv')# alternativelt save as csv, probably better
save(ScoresFDR_two.sided, file = "Table2_Meta40_FDR_default_TCGArem.rda") # same the unfiltered normalised expression matrix



# BH-FDR

## https://www.youtube.com/watch?v=DRRT4k7A2dM

b<- as.data.frame(ScoresFDR_two.sided)
class(b)
b$rawPbyZ <- 2*pnorm(-abs(b$zSco))
View(b)
b$FDR_BH <- p.adjust(b$rawPbyZ, "BH") # I dont think it is necessary to do fdr if we assuming to looking gene one by one
View(b)
write.csv(b, 'Table3_Meta40_FDR_BH_TCGArem.csv')# alternativelt save as csv, probably better



###################################################
### Pooled ES volcano
###################################################

GeneMeta_res<- as.data.frame(b)
GeneMeta_res$pooled_ES <- GeneMeta_res$MUvals
View(GeneMeta_res)

deg.data <- as.data.frame(GeneMeta_res[ , c(19,18)])
View(deg.data)
class(deg.data)
deg.data$logP <- -log10(deg.data$FDR_BH)

deg.data$Group = "Not significant"
deg.data$Group[which( (deg.data$FDR_BH < 0.05))] = "Significant"
#deg.data$Group[which( (deg.data$adj.P.Val < p.adjust.cutoff) & (deg.data$logFC < -log(fold.change.cutoff,2)) )] = "Down-regulated"
deg.data$Label = ""
deg.data$LabelSF = ""
deg.data <- deg.data[order(deg.data$FDR_BH), ]
up.genes <- head(rownames(deg.data)[which(deg.data$Group == "Significant")], 10)
#down.genes <- head(deg.data$Symbol[which(deg.data$Group == "Down-regulated")], 10)
deg.top10.genes <- (as.character(up.genes))
deg.data$Label[match(deg.top10.genes, rownames(deg.data))] <- deg.top10.genes
deg.data$LabelSF[match(deg.top10.genes, rownames(deg.data))] <- deg.top10.genes # NEED to run ' Extracting splicing factor info' below beforehand for 'degSF'
deg.data$Label_all <- rownames(deg.data)
View(deg.data)

library(ggplot2)
library(ggpubr)
library(cowplot)
p1 <- ggscatter(deg.data, x = "pooled_ES", y = "logP",
                color = "Group",
                shape = "Group",
                palette = c("#616161", "#000000"),#c("#BBBBBB", "#CC0000"),
                size = 2.5,
                label = deg.data$Label_all,
                #label = deg.data$LabelSF,
                font.label = 10,
                repel = T,
                xlab = "Pooled effect size (Hedges'g)",
                xlim = c(-0.3,1.3),
                ylab = "-log10(BH adjusted p-value)",
                ylim = c(0,38)) +
  theme_bw() +
  geom_hline(yintercept = -log(0.05,10), linetype="dashed") +
  ggtitle("Volcano plot of gene markers")
#geom_vline(xintercept = c(log(fold.change.cutoff,2),-log(fold.change.cutoff,2)), linetype="dashed")+
#ggtitle(paste("Volcano plot","(",colnames(fit2$contrasts),")"))
p1

## output plot
save_plot("volcano.pdf",p1,base_width = 12, base_height = 14)
save_plot("volcano.tiff",p1,base_width = 12, base_height = 14,
          compression = "lzw",dpi = 500)

#-----------------------------------------------------------------##
#                          Forest                                 ##
#-----------------------------------------------------------------##   
library(meta)


#https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/pooling-es.html
#https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/effects.html

lname <- load('Table1_Meta40_TCGArem.rda')
lname

DF <- as.data.frame(theScores)
class(DF)
View(DF)


sig_genes <- deg.data [deg.data$FDR_BH < 0.05, ]
View(sig_genes)


# LOOP for meta forest

meta.genes <- rownames(sig_genes)

for(gene in meta.genes) {
  
  GEO.study <- c('GSE21034','GSE70768', 'GSE94767', 'CIT')
  
  
  GSE21034_DF <- DF[rownames(DF) == gene, ]
  GSE21034_ES <- GSE21034_DF$Effect_Ex_1
  GSE21034_se <- sqrt(GSE21034_DF$EffectVar_Ex_1)
  #GSE117999_se <- sqrt(((10+10)/(10*10))+((GSE117999_ES)^2/(2*(10+10))))
  
  
  GSE70768_DF <- DF[rownames(DF) == gene, ]
  GSE70768_ES <- GSE70768_DF$Effect_Ex_2
  GSE70768_se <- sqrt(GSE70768_DF$EffectVar_Ex_2)
  #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
  
  GSE94767_DF <- DF[rownames(DF) == gene, ]
  GSE94767_ES <- GSE94767_DF$Effect_Ex_3
  GSE94767_se <- sqrt(GSE94767_DF$EffectVar_Ex_3)
  #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
  
  CIT_DF <- DF[rownames(DF) == gene, ]
  CIT_ES <- CIT_DF$Effect_Ex_4
  CIT_se <- sqrt(CIT_DF$EffectVar_Ex_4)
  #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
  
  
  ## Variable set-up
  
  TE <- c(GSE21034_ES , GSE70768_ES, GSE94767_ES, CIT_ES )
  seTE <- c(GSE21034_se, GSE70768_se, GSE94767_se, CIT_se)
  
  # Join the variables to create a data frame
  df2_gene <- data.frame(GEO.study,TE,seTE )
  df2_gene
  rownames(df2_gene) <- df2_gene$GEO.study
  
  data.frame(df2_gene, stringsAsFactors = TRUE)
  
  m.gen <- metagen(TE = TE,
                   seTE = seTE,
                   studlab = GEO.study,
                   data = df2_gene,
                   sm = "SMD",
                   #method.ci = 't',
                   method.smd = "Hedges",
                   fixed = FALSE,
                   random = TRUE,
                   method.tau = "DL",#!!!!!!!!!!! here we use DJ (DerSimonian and Laird estimator) because this is what GeneMeta uses!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   #method.tau = "REML",
                   hakn = FALSE, #if use false then metafor/GeneMeta/meta all matched!!!!!!!!! AND it uses Wald-type test
                   #hakn = TRUE,
                   title = gene)
  m.gen
  summary(m.gen)
  m.gen$pval.random
  
  pdf(file = paste0(gene,".meta_forest.pdf"), width = 14, height = 8)
  forest(m.gen, studlab = TRUE)
  dev.off()
  
  
}









##################################################################################################################################################################################
################################################### For four studies (GSE21034 removed) ##########################################################################################
##################################################################################################################################################################################

###################################################
### code chunk number 13: claculationAllinOne
###################################################
esets     <- list(eset_TCGA,
                  
                  eset_GSE70768,
                  eset_GSE94767,
                  eset_CIT)
# data.ER   <-pData(Nevins)[,"ER.status"]
# levels(data.ER) <- c(0,1)
# data.ER<- as.numeric(as.character(data.ER))
classes   <- list(condition_TCGA,
                  
                  condition_GSE70768,
                  condition_GSE94767,
                  condition_CIT)
#theScores <- zScores(esets,classes,useREM=FALSE,CombineExp=1:2)
theScores <- zScores(esets,classes,useREM=TRUE,CombineExp=1:4) # CombineExp can be set to combine wanted datasets
View(theScores)



###################################################
### code chunk number 14: show results
###################################################
View(theScores)
theScores[1:2,]

# output data
#write.table(theScores,"Table1_Meta_synovial.txt",row.names = T,col.names = NA,quote = F,sep="\t")
write.csv(theScores, 'Table1_Meta40_GSE21034rem.csv')# alternativelt save as csv, probably better
save(theScores, file = "Table1_Meta40_GSE21034rem.rda") # same the unfiltered normalised expression matrix



# FDR
set.seed(123)
ScoresFDR <- zScoreFDR(esets, classes, useREM=TRUE, nperm=50,CombineExp=1:4)
View(ScoresFDR$two.sided)
#View(ScoresFDR$pos)

ScoresFDR_two.sided <- ScoresFDR$two.sided
View(ScoresFDR_two.sided)
write.csv(ScoresFDR_two.sided, 'Table2_Meta40_FDR_default_GSE21034rem.csv')# alternativelt save as csv, probably better
save(ScoresFDR_two.sided, file = "Table2_Meta40_FDR_default_GSE21034rem.rda") # same the unfiltered normalised expression matrix



# BH-FDR

## https://www.youtube.com/watch?v=DRRT4k7A2dM

b<- as.data.frame(ScoresFDR_two.sided)
class(b)
b$rawPbyZ <- 2*pnorm(-abs(b$zSco))
View(b)
b$FDR_BH <- p.adjust(b$rawPbyZ, "BH") # I dont think it is necessary to do fdr if we assuming to looking gene one by one
View(b)
write.csv(b, 'Table3_Meta40_FDR_BH_GSE21034rem.csv')# alternativelt save as csv, probably better



###################################################
### Pooled ES volcano
###################################################

GeneMeta_res<- as.data.frame(b)
GeneMeta_res$pooled_ES <- GeneMeta_res$MUvals
View(GeneMeta_res)

deg.data <- as.data.frame(GeneMeta_res[ , c(19,18)])
View(deg.data)
class(deg.data)
deg.data$logP <- -log10(deg.data$FDR_BH)

deg.data$Group = "Not significant"
deg.data$Group[which( (deg.data$FDR_BH < 0.05))] = "Significant"
#deg.data$Group[which( (deg.data$adj.P.Val < p.adjust.cutoff) & (deg.data$logFC < -log(fold.change.cutoff,2)) )] = "Down-regulated"
deg.data$Label = ""
deg.data$LabelSF = ""
deg.data <- deg.data[order(deg.data$FDR_BH), ]
up.genes <- head(rownames(deg.data)[which(deg.data$Group == "Significant")], 10)
#down.genes <- head(deg.data$Symbol[which(deg.data$Group == "Down-regulated")], 10)
deg.top10.genes <- (as.character(up.genes))
deg.data$Label[match(deg.top10.genes, rownames(deg.data))] <- deg.top10.genes
deg.data$LabelSF[match(deg.top10.genes, rownames(deg.data))] <- deg.top10.genes # NEED to run ' Extracting splicing factor info' below beforehand for 'degSF'
deg.data$Label_all <- rownames(deg.data)
View(deg.data)

library(ggplot2)
library(ggpubr)
p1 <- ggscatter(deg.data, x = "pooled_ES", y = "logP",
                color = "Group",
                shape = "Group",
                palette = c("#616161", "#000000"),#c("#BBBBBB", "#CC0000"),
                size = 2.5,
                label = deg.data$Label_all,
                #label = deg.data$LabelSF,
                font.label = 10,
                repel = T,
                xlab = "Pooled effect size (Hedges'g)",
                xlim = c(-0.5,1.4),
                ylab = "-log10(BH adjusted p-value)",
                ylim = c(0,55)) +
  theme_bw() +
  geom_hline(yintercept = -log(0.05,10), linetype="dashed") +
  ggtitle("Volcano plot of gene markers")
#geom_vline(xintercept = c(log(fold.change.cutoff,2),-log(fold.change.cutoff,2)), linetype="dashed")+
#ggtitle(paste("Volcano plot","(",colnames(fit2$contrasts),")"))
p1
library(cowplot)
## output plot
save_plot("volcano.pdf",p1,base_width = 12, base_height = 14)
save_plot("volcano.tiff",p1,base_width = 12, base_height = 14,
          compression = "lzw",dpi = 500)

#-----------------------------------------------------------------##
#                          Forest                                 ##
#-----------------------------------------------------------------##   
library(meta)


#https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/pooling-es.html
#https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/effects.html

lname <- load('Table1_Meta40_GSE21034rem.rda')
lname

DF <- as.data.frame(theScores)
class(DF)
View(DF)


sig_genes <- deg.data [deg.data$FDR_BH < 0.05, ]
View(sig_genes)


# LOOP for meta forest

meta.genes <- rownames(sig_genes)

for(gene in meta.genes) {
  
  GEO.study <- c('TCGA','GSE70768', 'GSE94767', 'CIT')
  
  TCGA_DF <- DF[rownames(DF) == gene, ]
  TCGA_ES <- TCGA_DF$Effect_Ex_1
  TCGA_se <- sqrt(TCGA_DF$EffectVar_Ex_1)
  #GSE57218_se <- sqrt(((5+5)/(5*5))+((GSE57218_ES)^2/(2*(5+5))))
  

  
  
  GSE70768_DF <- DF[rownames(DF) == gene, ]
  GSE70768_ES <- GSE70768_DF$Effect_Ex_2
  GSE70768_se <- sqrt(GSE70768_DF$EffectVar_Ex_2)
  #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
  
  GSE94767_DF <- DF[rownames(DF) == gene, ]
  GSE94767_ES <- GSE94767_DF$Effect_Ex_3
  GSE94767_se <- sqrt(GSE94767_DF$EffectVar_Ex_3)
  #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
  
  CIT_DF <- DF[rownames(DF) == gene, ]
  CIT_ES <- CIT_DF$Effect_Ex_4
  CIT_se <- sqrt(CIT_DF$EffectVar_Ex_4)
  #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
  
  
  ## Variable set-up
  
  TE <- c(TCGA_ES,  GSE70768_ES, GSE94767_ES, CIT_ES )
  seTE <- c(TCGA_se, GSE70768_se, GSE94767_se, CIT_se)
  
  # Join the variables to create a data frame
  df2_gene <- data.frame(GEO.study,TE,seTE )
  df2_gene
  rownames(df2_gene) <- df2_gene$GEO.study
  
  data.frame(df2_gene, stringsAsFactors = TRUE)
  
  m.gen <- metagen(TE = TE,
                   seTE = seTE,
                   studlab = GEO.study,
                   data = df2_gene,
                   sm = "SMD",
                   #method.ci = 't',
                   method.smd = "Hedges",
                   fixed = FALSE,
                   random = TRUE,
                   method.tau = "DL",#!!!!!!!!!!! here we use DJ (DerSimonian and Laird estimator) because this is what GeneMeta uses!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   #method.tau = "REML",
                   hakn = FALSE, #if use false then metafor/GeneMeta/meta all matched!!!!!!!!! AND it uses Wald-type test
                   #hakn = TRUE,
                   title = gene)
  m.gen
  summary(m.gen)
  m.gen$pval.random
  
  pdf(file = paste0(gene,".meta_forest.pdf"), width = 14, height = 8)
  forest(m.gen, studlab = TRUE)
  dev.off()
  
  
}












##################################################################################################################################################################################
################################################################## For four studies GSE70768 removed #############################################################################
##################################################################################################################################################################################

###################################################
### code chunk number 13: claculationAllinOne
###################################################
esets     <- list(eset_TCGA,
                  eset_GSE21034,
                  
                  eset_GSE94767,
                  eset_CIT)
# data.ER   <-pData(Nevins)[,"ER.status"]
# levels(data.ER) <- c(0,1)
# data.ER<- as.numeric(as.character(data.ER))
classes   <- list(condition_TCGA,
                  condition_GSE21034,
                  
                  condition_GSE94767,
                  condition_CIT)
#theScores <- zScores(esets,classes,useREM=FALSE,CombineExp=1:2)
theScores <- zScores(esets,classes,useREM=TRUE,CombineExp=1:4) # CombineExp can be set to combine wanted datasets
View(theScores)



###################################################
### code chunk number 14: show results
###################################################
View(theScores)
theScores[1:2,]

# output data
#write.table(theScores,"Table1_Meta_synovial.txt",row.names = T,col.names = NA,quote = F,sep="\t")
write.csv(theScores, 'Table1_Meta40_GSE70768rem.csv')# alternativelt save as csv, probably better
save(theScores, file = "Table1_Meta40_GSE70768rem.rda") # same the unfiltered normalised expression matrix



# FDR
set.seed(123)
ScoresFDR <- zScoreFDR(esets, classes, useREM=TRUE, nperm=50,CombineExp=1:4)
View(ScoresFDR$two.sided)
#View(ScoresFDR$pos)

ScoresFDR_two.sided <- ScoresFDR$two.sided
View(ScoresFDR_two.sided)
write.csv(ScoresFDR_two.sided, 'Table2_Meta40_FDR_default_GSE70768rem.csv')# alternativelt save as csv, probably better
save(ScoresFDR_two.sided, file = "Table2_Meta40_FDR_default_GSE70768rem.rda") # same the unfiltered normalised expression matrix



# BH-FDR

## https://www.youtube.com/watch?v=DRRT4k7A2dM

b<- as.data.frame(ScoresFDR_two.sided)
class(b)
b$rawPbyZ <- 2*pnorm(-abs(b$zSco))
View(b)
b$FDR_BH <- p.adjust(b$rawPbyZ, "BH") # I dont think it is necessary to do fdr if we assuming to looking gene one by one
View(b)
write.csv(b, 'Table3_Meta40_FDR_BH_GSE70768rem.csv')# alternativelt save as csv, probably better



###################################################
### Pooled ES volcano
###################################################

GeneMeta_res<- as.data.frame(b)
GeneMeta_res$pooled_ES <- GeneMeta_res$MUvals
View(GeneMeta_res)

deg.data <- as.data.frame(GeneMeta_res[ , c(19,18)])
View(deg.data)
class(deg.data)
deg.data$logP <- -log10(deg.data$FDR_BH)

deg.data$Group = "Not significant"
deg.data$Group[which( (deg.data$FDR_BH < 0.05))] = "Significant"
#deg.data$Group[which( (deg.data$adj.P.Val < p.adjust.cutoff) & (deg.data$logFC < -log(fold.change.cutoff,2)) )] = "Down-regulated"
deg.data$Label = ""
deg.data$LabelSF = ""
deg.data <- deg.data[order(deg.data$FDR_BH), ]
up.genes <- head(rownames(deg.data)[which(deg.data$Group == "Significant")], 10)
#down.genes <- head(deg.data$Symbol[which(deg.data$Group == "Down-regulated")], 10)
deg.top10.genes <- (as.character(up.genes))
deg.data$Label[match(deg.top10.genes, rownames(deg.data))] <- deg.top10.genes
deg.data$LabelSF[match(deg.top10.genes, rownames(deg.data))] <- deg.top10.genes # NEED to run ' Extracting splicing factor info' below beforehand for 'degSF'
deg.data$Label_all <- rownames(deg.data)
View(deg.data)

library(ggplot2)
library(ggpubr)
p1 <- ggscatter(deg.data, x = "pooled_ES", y = "logP",
                color = "Group",
                shape = "Group",
                palette = c("#616161", "#000000"),#c("#BBBBBB", "#CC0000"),
                size = 2.5,
                label = deg.data$Label_all,
                #label = deg.data$LabelSF,
                font.label = 10,
                repel = T,
                xlab = "Pooled effect size (Hedges'g)",
                xlim = c(-0.4,1.3),
                ylab = "-log10(BH adjusted p-value)",
                ylim = c(0,40)) +
  theme_bw() +
  geom_hline(yintercept = -log(0.05,10), linetype="dashed") +
  ggtitle("Volcano plot of gene markers")
#geom_vline(xintercept = c(log(fold.change.cutoff,2),-log(fold.change.cutoff,2)), linetype="dashed")+
#ggtitle(paste("Volcano plot","(",colnames(fit2$contrasts),")"))
p1
library(cowplot)
## output plot
save_plot("volcano.pdf",p1,base_width = 12, base_height = 14)
save_plot("volcano.tiff",p1,base_width = 12, base_height = 14,
          compression = "lzw",dpi = 500)

#-----------------------------------------------------------------##
#                          Forest                                 ##
#-----------------------------------------------------------------##   
library(meta)


#https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/pooling-es.html
#https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/effects.html

lname <- load('Table1_Meta40_GSE70768rem.rda')
lname

DF <- as.data.frame(theScores)
class(DF)
View(DF)


sig_genes <- deg.data [deg.data$FDR_BH < 0.05, ]
View(sig_genes)


# LOOP for meta forest

meta.genes <- rownames(sig_genes)

for(gene in meta.genes) {
  
  GEO.study <- c('TCGA','GSE21034', 'GSE94767', 'CIT')
  
  TCGA_DF <- DF[rownames(DF) == gene, ]
  TCGA_ES <- TCGA_DF$Effect_Ex_1
  TCGA_se <- sqrt(TCGA_DF$EffectVar_Ex_1)
  #GSE57218_se <- sqrt(((5+5)/(5*5))+((GSE57218_ES)^2/(2*(5+5))))
  
  GSE21034_DF <- DF[rownames(DF) == gene, ]
  GSE21034_ES <- GSE21034_DF$Effect_Ex_2
  GSE21034_se <- sqrt(GSE21034_DF$EffectVar_Ex_2)
  #GSE117999_se <- sqrt(((10+10)/(10*10))+((GSE117999_ES)^2/(2*(10+10))))
  
  
  
  GSE94767_DF <- DF[rownames(DF) == gene, ]
  GSE94767_ES <- GSE94767_DF$Effect_Ex_3
  GSE94767_se <- sqrt(GSE94767_DF$EffectVar_Ex_3)
  #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
  
  CIT_DF <- DF[rownames(DF) == gene, ]
  CIT_ES <- CIT_DF$Effect_Ex_4
  CIT_se <- sqrt(CIT_DF$EffectVar_Ex_4)
  #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
  
  
  ## Variable set-up
  
  TE <- c(TCGA_ES, GSE21034_ES , GSE94767_ES, CIT_ES )
  seTE <- c(TCGA_se,GSE21034_se, GSE94767_se, CIT_se)
  
  # Join the variables to create a data frame
  df2_gene <- data.frame(GEO.study,TE,seTE )
  df2_gene
  rownames(df2_gene) <- df2_gene$GEO.study
  
  data.frame(df2_gene, stringsAsFactors = TRUE)
  
  m.gen <- metagen(TE = TE,
                   seTE = seTE,
                   studlab = GEO.study,
                   data = df2_gene,
                   sm = "SMD",
                   #method.ci = 't',
                   method.smd = "Hedges",
                   fixed = FALSE,
                   random = TRUE,
                   method.tau = "DL",#!!!!!!!!!!! here we use DJ (DerSimonian and Laird estimator) because this is what GeneMeta uses!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   #method.tau = "REML",
                   hakn = FALSE, #if use false then metafor/GeneMeta/meta all matched!!!!!!!!! AND it uses Wald-type test
                   #hakn = TRUE,
                   title = gene)
  m.gen
  summary(m.gen)
  m.gen$pval.random
  
  pdf(file = paste0(gene,".meta_forest.pdf"), width = 14, height = 8)
  forest(m.gen, studlab = TRUE)
  dev.off()
  
  
}







##################################################################################################################################################################################
################################################################## For four studies GSE94767 removed ############################################################################
##################################################################################################################################################################################

###################################################
### code chunk number 13: claculationAllinOne
###################################################
esets     <- list(eset_TCGA,
                  eset_GSE21034,
                  eset_GSE70768,
                
                  eset_CIT)
# data.ER   <-pData(Nevins)[,"ER.status"]
# levels(data.ER) <- c(0,1)
# data.ER<- as.numeric(as.character(data.ER))
classes   <- list(condition_TCGA,
                  condition_GSE21034,
                  condition_GSE70768,
                  
                  condition_CIT)
#theScores <- zScores(esets,classes,useREM=FALSE,CombineExp=1:2)
theScores <- zScores(esets,classes,useREM=TRUE,CombineExp=1:4) # CombineExp can be set to combine wanted datasets
View(theScores)



###################################################
### code chunk number 14: show results
###################################################
View(theScores)
theScores[1:2,]

# output data
#write.table(theScores,"Table1_Meta_synovial.txt",row.names = T,col.names = NA,quote = F,sep="\t")
write.csv(theScores, 'Table1_Meta40_GSE94767rem.csv')# alternativelt save as csv, probably better
save(theScores, file = "Table1_Meta40_GSE94767rem.rda") # same the unfiltered normalised expression matrix



# FDR
set.seed(123)
ScoresFDR <- zScoreFDR(esets, classes, useREM=TRUE, nperm=50,CombineExp=1:4)
View(ScoresFDR$two.sided)
#View(ScoresFDR$pos)

ScoresFDR_two.sided <- ScoresFDR$two.sided
View(ScoresFDR_two.sided)
write.csv(ScoresFDR_two.sided, 'Table2_Meta40_FDR_default_GSE94767rem.csv')# alternativelt save as csv, probably better
save(ScoresFDR_two.sided, file = "Table2_Meta40_FDR_default_GSE94767rem.rda") # same the unfiltered normalised expression matrix



# BH-FDR

## https://www.youtube.com/watch?v=DRRT4k7A2dM

b<- as.data.frame(ScoresFDR_two.sided)
class(b)
b$rawPbyZ <- 2*pnorm(-abs(b$zSco))
View(b)
b$FDR_BH <- p.adjust(b$rawPbyZ, "BH") # I dont think it is necessary to do fdr if we assuming to looking gene one by one
View(b)
write.csv(b, 'Table3_Meta40_FDR_BH_GSE94767rem.csv')# alternativelt save as csv, probably better



###################################################
### Pooled ES volcano
###################################################

GeneMeta_res<- as.data.frame(b)
GeneMeta_res$pooled_ES <- GeneMeta_res$MUvals
View(GeneMeta_res)

deg.data <- as.data.frame(GeneMeta_res[ , c(19,18)])
View(deg.data)
class(deg.data)
deg.data$logP <- -log10(deg.data$FDR_BH)

deg.data$Group = "Not significant"
deg.data$Group[which( (deg.data$FDR_BH < 0.05))] = "Significant"
#deg.data$Group[which( (deg.data$adj.P.Val < p.adjust.cutoff) & (deg.data$logFC < -log(fold.change.cutoff,2)) )] = "Down-regulated"
deg.data$Label = ""
deg.data$LabelSF = ""
deg.data <- deg.data[order(deg.data$FDR_BH), ]
up.genes <- head(rownames(deg.data)[which(deg.data$Group == "Significant")], 10)
#down.genes <- head(deg.data$Symbol[which(deg.data$Group == "Down-regulated")], 10)
deg.top10.genes <- (as.character(up.genes))
deg.data$Label[match(deg.top10.genes, rownames(deg.data))] <- deg.top10.genes
deg.data$LabelSF[match(deg.top10.genes, rownames(deg.data))] <- deg.top10.genes # NEED to run ' Extracting splicing factor info' below beforehand for 'degSF'
deg.data$Label_all <- rownames(deg.data)
View(deg.data)

library(ggplot2)
library(ggpubr)
p1 <- ggscatter(deg.data, x = "pooled_ES", y = "logP",
                color = "Group",
                shape = "Group",
                palette = c("#616161", "#000000"),#c("#BBBBBB", "#CC0000"),
                size = 2.5,
                label = deg.data$Label_all,
                #label = deg.data$LabelSF,
                font.label = 10,
                repel = T,
                xlab = "Pooled effect size (Hedges'g)",
                xlim = c(-0.4,1.3),
                ylab = "-log10(BH adjusted p-value)",
                ylim = c(0,50)) +
  theme_bw() +
  geom_hline(yintercept = -log(0.05,10), linetype="dashed") +
  ggtitle("Volcano plot of gene markers")
#geom_vline(xintercept = c(log(fold.change.cutoff,2),-log(fold.change.cutoff,2)), linetype="dashed")+
#ggtitle(paste("Volcano plot","(",colnames(fit2$contrasts),")"))
p1
library(cowplot)
## output plot
save_plot("volcano.pdf",p1,base_width = 12, base_height = 14)
save_plot("volcano.tiff",p1,base_width = 12, base_height = 14,
          compression = "lzw",dpi = 500)

#-----------------------------------------------------------------##
#                          Forest                                 ##
#-----------------------------------------------------------------##   
library(meta)


#https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/pooling-es.html
#https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/effects.html

lname <- load('Table1_Meta40_GSE94767rem.rda')
lname

DF <- as.data.frame(theScores)
class(DF)
View(DF)


sig_genes <- deg.data [deg.data$FDR_BH < 0.05, ]
View(sig_genes)


# LOOP for meta forest

meta.genes <- rownames(sig_genes)

for(gene in meta.genes) {
  
  GEO.study <- c('TCGA','GSE21034','GSE70768', 'CIT')
  
  TCGA_DF <- DF[rownames(DF) == gene, ]
  TCGA_ES <- TCGA_DF$Effect_Ex_1
  TCGA_se <- sqrt(TCGA_DF$EffectVar_Ex_1)
  #GSE57218_se <- sqrt(((5+5)/(5*5))+((GSE57218_ES)^2/(2*(5+5))))
  
  GSE21034_DF <- DF[rownames(DF) == gene, ]
  GSE21034_ES <- GSE21034_DF$Effect_Ex_2
  GSE21034_se <- sqrt(GSE21034_DF$EffectVar_Ex_2)
  #GSE117999_se <- sqrt(((10+10)/(10*10))+((GSE117999_ES)^2/(2*(10+10))))
  
  
  GSE70768_DF <- DF[rownames(DF) == gene, ]
  GSE70768_ES <- GSE70768_DF$Effect_Ex_3
  GSE70768_se <- sqrt(GSE70768_DF$EffectVar_Ex_3)
  #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
  

  CIT_DF <- DF[rownames(DF) == gene, ]
  CIT_ES <- CIT_DF$Effect_Ex_4
  CIT_se <- sqrt(CIT_DF$EffectVar_Ex_4)
  #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
  
  
  ## Variable set-up
  
  TE <- c(TCGA_ES, GSE21034_ES , GSE70768_ES,  CIT_ES )
  seTE <- c(TCGA_se,GSE21034_se, GSE70768_se, CIT_se)
  
  # Join the variables to create a data frame
  df2_gene <- data.frame(GEO.study,TE,seTE )
  df2_gene
  rownames(df2_gene) <- df2_gene$GEO.study
  
  data.frame(df2_gene, stringsAsFactors = TRUE)
  
  m.gen <- metagen(TE = TE,
                   seTE = seTE,
                   studlab = GEO.study,
                   data = df2_gene,
                   sm = "SMD",
                   #method.ci = 't',
                   method.smd = "Hedges",
                   fixed = FALSE,
                   random = TRUE,
                   method.tau = "DL",#!!!!!!!!!!! here we use DJ (DerSimonian and Laird estimator) because this is what GeneMeta uses!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   #method.tau = "REML",
                   hakn = FALSE, #if use false then metafor/GeneMeta/meta all matched!!!!!!!!! AND it uses Wald-type test
                   #hakn = TRUE,
                   title = gene)
  m.gen
  summary(m.gen)
  m.gen$pval.random
  
  pdf(file = paste0(gene,".meta_forest.pdf"), width = 14, height = 8)
  forest(m.gen, studlab = TRUE)
  dev.off()
  
  
}









##################################################################################################################################################################################
################################################################## For four studies CIT removed ##################################################################################
##################################################################################################################################################################################

###################################################
### code chunk number 13: claculationAllinOne
###################################################
esets     <- list(eset_TCGA,
                  eset_GSE21034,
                  eset_GSE70768,
                  eset_GSE94767)
# data.ER   <-pData(Nevins)[,"ER.status"]
# levels(data.ER) <- c(0,1)
# data.ER<- as.numeric(as.character(data.ER))
classes   <- list(condition_TCGA,
                  condition_GSE21034,
                  condition_GSE70768,
                  condition_GSE94767)
#theScores <- zScores(esets,classes,useREM=FALSE,CombineExp=1:2)
theScores <- zScores(esets,classes,useREM=TRUE,CombineExp=1:4) # CombineExp can be set to combine wanted datasets
View(theScores)



###################################################
### code chunk number 14: show results
###################################################
View(theScores)
theScores[1:2,]

# output data
#write.table(theScores,"Table1_Meta_synovial.txt",row.names = T,col.names = NA,quote = F,sep="\t")
write.csv(theScores, 'Table1_Meta40_CITrem.csv')# alternativelt save as csv, probably better
save(theScores, file = "Table1_Meta40_CITrem.rda") # same the unfiltered normalised expression matrix



# FDR
set.seed(123)
ScoresFDR <- zScoreFDR(esets, classes, useREM=TRUE, nperm=50,CombineExp=1:4)
View(ScoresFDR$two.sided)
#View(ScoresFDR$pos)

ScoresFDR_two.sided <- ScoresFDR$two.sided
View(ScoresFDR_two.sided)
write.csv(ScoresFDR_two.sided, 'Table2_Meta40_FDR_default_CITrem.csv')# alternativelt save as csv, probably better
save(ScoresFDR_two.sided, file = "Table2_Meta40_FDR_default_CITrem.rda") # same the unfiltered normalised expression matrix



# BH-FDR

## https://www.youtube.com/watch?v=DRRT4k7A2dM

b<- as.data.frame(ScoresFDR_two.sided)
class(b)
b$rawPbyZ <- 2*pnorm(-abs(b$zSco))
View(b)
b$FDR_BH <- p.adjust(b$rawPbyZ, "BH") # I dont think it is necessary to do fdr if we assuming to looking gene one by one
View(b)
write.csv(b, 'Table3_Meta40_FDR_BH_CITrem.csv')# alternativelt save as csv, probably better



###################################################
### Pooled ES volcano
###################################################

GeneMeta_res<- as.data.frame(b)
GeneMeta_res$pooled_ES <- GeneMeta_res$MUvals
View(GeneMeta_res)

deg.data <- as.data.frame(GeneMeta_res[ , c(19,18)])
View(deg.data)
class(deg.data)
deg.data$logP <- -log10(deg.data$FDR_BH)

deg.data$Group = "Not significant"
deg.data$Group[which( (deg.data$FDR_BH < 0.05))] = "Significant"
#deg.data$Group[which( (deg.data$adj.P.Val < p.adjust.cutoff) & (deg.data$logFC < -log(fold.change.cutoff,2)) )] = "Down-regulated"
deg.data$Label = ""
deg.data$LabelSF = ""
deg.data <- deg.data[order(deg.data$FDR_BH), ]
up.genes <- head(rownames(deg.data)[which(deg.data$Group == "Significant")], 10)
#down.genes <- head(deg.data$Symbol[which(deg.data$Group == "Down-regulated")], 10)
deg.top10.genes <- (as.character(up.genes))
deg.data$Label[match(deg.top10.genes, rownames(deg.data))] <- deg.top10.genes
deg.data$LabelSF[match(deg.top10.genes, rownames(deg.data))] <- deg.top10.genes # NEED to run ' Extracting splicing factor info' below beforehand for 'degSF'
deg.data$Label_all <- rownames(deg.data)
View(deg.data)

library(ggplot2)
library(ggpubr)
p1 <- ggscatter(deg.data, x = "pooled_ES", y = "logP",
                color = "Group",
                shape = "Group",
                palette = c("#616161", "#000000"),#c("#BBBBBB", "#CC0000"),
                size = 2.5,
                label = deg.data$Label_all,
                #label = deg.data$LabelSF,
                font.label = 10,
                repel = T,
                xlab = "Pooled effect size (Hedges'g)",
                xlim = c(-0.5,1.4),
                ylab = "-log10(BH adjusted p-value)",
                ylim = c(0,50)) +
  theme_bw() +
  geom_hline(yintercept = -log(0.05,10), linetype="dashed") +
  ggtitle("Volcano plot of gene markers")
#geom_vline(xintercept = c(log(fold.change.cutoff,2),-log(fold.change.cutoff,2)), linetype="dashed")+
#ggtitle(paste("Volcano plot","(",colnames(fit2$contrasts),")"))
p1
library(cowplot)
## output plot
save_plot("volcano.pdf",p1,base_width = 12, base_height = 14)
save_plot("volcano.tiff",p1,base_width = 12, base_height = 14,
          compression = "lzw",dpi = 500)

#-----------------------------------------------------------------##
#                          Forest                                 ##
#-----------------------------------------------------------------##   
library(meta)


#https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/pooling-es.html
#https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/effects.html

lname <- load('Table1_Meta40_CITrem.rda')
lname

DF <- as.data.frame(theScores)
class(DF)
View(DF)


sig_genes <- deg.data [deg.data$FDR_BH < 0.05, ]
View(sig_genes)


# LOOP for meta forest

meta.genes <- rownames(sig_genes)

for(gene in meta.genes) {
  
  GEO.study <- c('TCGA','GSE21034','GSE70768', 'GSE94767')
  
  TCGA_DF <- DF[rownames(DF) == gene, ]
  TCGA_ES <- TCGA_DF$Effect_Ex_1
  TCGA_se <- sqrt(TCGA_DF$EffectVar_Ex_1)
  #GSE57218_se <- sqrt(((5+5)/(5*5))+((GSE57218_ES)^2/(2*(5+5))))
  
  GSE21034_DF <- DF[rownames(DF) == gene, ]
  GSE21034_ES <- GSE21034_DF$Effect_Ex_2
  GSE21034_se <- sqrt(GSE21034_DF$EffectVar_Ex_2)
  #GSE117999_se <- sqrt(((10+10)/(10*10))+((GSE117999_ES)^2/(2*(10+10))))
  
  
  GSE70768_DF <- DF[rownames(DF) == gene, ]
  GSE70768_ES <- GSE70768_DF$Effect_Ex_3
  GSE70768_se <- sqrt(GSE70768_DF$EffectVar_Ex_3)
  #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
  
  GSE94767_DF <- DF[rownames(DF) == gene, ]
  GSE94767_ES <- GSE94767_DF$Effect_Ex_4
  GSE94767_se <- sqrt(GSE94767_DF$EffectVar_Ex_4)
  #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
  

  
  
  ## Variable set-up
  
  TE <- c(TCGA_ES, GSE21034_ES , GSE70768_ES, GSE94767_ES )
  seTE <- c(TCGA_se,GSE21034_se, GSE70768_se, GSE94767_se)
  
  # Join the variables to create a data frame
  df2_gene <- data.frame(GEO.study,TE,seTE )
  df2_gene
  rownames(df2_gene) <- df2_gene$GEO.study
  
  data.frame(df2_gene, stringsAsFactors = TRUE)
  
  m.gen <- metagen(TE = TE,
                   seTE = seTE,
                   studlab = GEO.study,
                   data = df2_gene,
                   sm = "SMD",
                   #method.ci = 't',
                   method.smd = "Hedges",
                   fixed = FALSE,
                   random = TRUE,
                   method.tau = "DL",#!!!!!!!!!!! here we use DJ (DerSimonian and Laird estimator) because this is what GeneMeta uses!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   #method.tau = "REML",
                   hakn = FALSE, #if use false then metafor/GeneMeta/meta all matched!!!!!!!!! AND it uses Wald-type test
                   #hakn = TRUE,
                   title = gene)
  m.gen
  summary(m.gen)
  m.gen$pval.random
  
  pdf(file = paste0(gene,".meta_forest.pdf"), width = 14, height = 8)
  forest(m.gen, studlab = TRUE)
  dev.off()
  
  
}



# ## FBP1
# GEO.study <- c('TCGA','GSE21034','GSE70768', 'GSE94767', 'CIT')
# 
# TCGA_DF <- DF[rownames(DF) == 'FBP1', ]
# TCGA_ES <- TCGA_DF$Effect_Ex_1
# TCGA_se <- sqrt(TCGA_DF$EffectVar_Ex_1)
# #GSE57218_se <- sqrt(((5+5)/(5*5))+((GSE57218_ES)^2/(2*(5+5))))
# 
# GSE21034_DF <- DF[rownames(DF) == 'FBP1', ]
# GSE21034_ES <- GSE21034_DF$Effect_Ex_2
# GSE21034_se <- sqrt(GSE21034_DF$EffectVar_Ex_2)
# #GSE117999_se <- sqrt(((10+10)/(10*10))+((GSE117999_ES)^2/(2*(10+10))))
# 
# 
# GSE70768_DF <- DF[rownames(DF) == 'FBP1', ]
# GSE70768_ES <- GSE70768_DF$Effect_Ex_3
# GSE70768_se <- sqrt(GSE70768_DF$EffectVar_Ex_3)
# #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
# 
# GSE94767_DF <- DF[rownames(DF) == 'FBP1', ]
# GSE94767_ES <- GSE94767_DF$Effect_Ex_4
# GSE94767_se <- sqrt(GSE94767_DF$EffectVar_Ex_4)
# #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
# 
# CIT_DF <- DF[rownames(DF) == 'FBP1', ]
# CIT_ES <- CIT_DF$Effect_Ex_5
# CIT_se <- sqrt(CIT_DF$EffectVar_Ex_5)
# #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
# 
# 
# ## Variable set-up
# 
# TE <- c(TCGA_ES, GSE21034_ES , GSE70768_ES, GSE94767_ES, CIT_ES )
# seTE <- c(TCGA_se,GSE21034_se, GSE70768_se, GSE94767_se, CIT_se)
# 
# # Join the variables to create a data frame
# df2_FBP1 <- data.frame(GEO.study,TE,seTE )
# df2_FBP1
# rownames(df2_FBP1) <- df2_FBP1$GEO.study
# 
# data.frame(df2_FBP1, stringsAsFactors = TRUE)
# class(df2_FBP1)
# str(df2_FBP1)
# View(df2_FBP1)
# 
# 
# m.gen <- metagen(TE = TE,
#                  seTE = seTE,
#                  studlab = GEO.study,
#                  data = df2_FBP1,
#                  sm = "SMD",
#                  #method.ci = 't',
#                  method.smd = "Hedges",
#                  fixed = FALSE,
#                  random = TRUE,
#                  method.tau = "DL",#!!!!!!!!!!! here we use DJ (DerSimonian and Laird estimator) because this is what GeneMeta uses!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#                  #method.tau = "REML",
#                  hakn = FALSE, #if use false then metafor/GeneMeta/meta all matched!!!!!!!!! AND it uses Wald-type test
#                  #hakn = TRUE,
#                  title = "FBP1")
# m.gen
# summary(m.gen)
# m.gen$pval.random
# 
# pdf(file = "meta_forest_FBP1.pdf", width = 14, height = 8)
# forest(m.gen, studlab = TRUE)
# dev.off()
# 
# 
# 
# 
# ## MYC
# GEO.study <- c('TCGA','GSE21034','GSE70768', 'GSE94767', 'CIT')
# 
# TCGA_DF <- DF[rownames(DF) == 'MYC', ]
# TCGA_ES <- TCGA_DF$Effect_Ex_1
# TCGA_se <- sqrt(TCGA_DF$EffectVar_Ex_1)
# #GSE57218_se <- sqrt(((5+5)/(5*5))+((GSE57218_ES)^2/(2*(5+5))))
# 
# GSE21034_DF <- DF[rownames(DF) == 'MYC', ]
# GSE21034_ES <- GSE21034_DF$Effect_Ex_2
# GSE21034_se <- sqrt(GSE21034_DF$EffectVar_Ex_2)
# #GSE117999_se <- sqrt(((10+10)/(10*10))+((GSE117999_ES)^2/(2*(10+10))))
# 
# 
# GSE70768_DF <- DF[rownames(DF) == 'MYC', ]
# GSE70768_ES <- GSE70768_DF$Effect_Ex_3
# GSE70768_se <- sqrt(GSE70768_DF$EffectVar_Ex_3)
# #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
# 
# GSE94767_DF <- DF[rownames(DF) == 'MYC', ]
# GSE94767_ES <- GSE94767_DF$Effect_Ex_4
# GSE94767_se <- sqrt(GSE94767_DF$EffectVar_Ex_4)
# #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
# 
# CIT_DF <- DF[rownames(DF) == 'MYC', ]
# CIT_ES <- CIT_DF$Effect_Ex_5
# CIT_se <- sqrt(CIT_DF$EffectVar_Ex_5)
# #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
# 
# 
# ## Variable set-up
# 
# TE <- c(TCGA_ES, GSE21034_ES , GSE70768_ES, GSE94767_ES, CIT_ES )
# seTE <- c(TCGA_se,GSE21034_se, GSE70768_se, GSE94767_se, CIT_se)
# 
# # Join the variables to create a data frame
# df2_MYC <- data.frame(GEO.study,TE,seTE )
# df2_MYC
# rownames(df2_MYC) <- df2_MYC$GEO.study
# 
# data.frame(df2_MYC, stringsAsFactors = TRUE)
# class(df2_MYC)
# str(df2_MYC)
# View(df2_MYC)
# 
# 
# m.gen <- metagen(TE = TE,
#                  seTE = seTE,
#                  studlab = GEO.study,
#                  data = df2_MYC,
#                  sm = "SMD",
#                  #method.ci = 't',
#                  method.smd = "Hedges",
#                  fixed = FALSE,
#                  random = TRUE,
#                  method.tau = "DL",#!!!!!!!!!!! here we use DJ (DerSimonian and Laird estimator) because this is what GeneMeta uses!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#                  #method.tau = "REML",
#                  hakn = FALSE, #if use false then metafor/GeneMeta/meta all matched!!!!!!!!! AND it uses Wald-type test
#                  #hakn = TRUE,
#                  title = "MYC")
# m.gen
# summary(m.gen)
# m.gen$pval.random
# 
# pdf(file = "meta_forest_MYC.pdf", width = 14, height = 8)
# forest(m.gen, studlab = TRUE)
# dev.off()
# 
# 
# 
# 
# ## LAMTOR2
# GEO.study <- c('TCGA','GSE21034','GSE70768', 'GSE94767', 'CIT')
# 
# TCGA_DF <- DF[rownames(DF) == 'LAMTOR2', ]
# TCGA_ES <- TCGA_DF$Effect_Ex_1
# TCGA_se <- sqrt(TCGA_DF$EffectVar_Ex_1)
# #GSE57218_se <- sqrt(((5+5)/(5*5))+((GSE57218_ES)^2/(2*(5+5))))
# 
# GSE21034_DF <- DF[rownames(DF) == 'LAMTOR2', ]
# GSE21034_ES <- GSE21034_DF$Effect_Ex_2
# GSE21034_se <- sqrt(GSE21034_DF$EffectVar_Ex_2)
# #GSE117999_se <- sqrt(((10+10)/(10*10))+((GSE117999_ES)^2/(2*(10+10))))
# 
# 
# GSE70768_DF <- DF[rownames(DF) == 'LAMTOR2', ]
# GSE70768_ES <- GSE70768_DF$Effect_Ex_3
# GSE70768_se <- sqrt(GSE70768_DF$EffectVar_Ex_3)
# #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
# 
# GSE94767_DF <- DF[rownames(DF) == 'LAMTOR2', ]
# GSE94767_ES <- GSE94767_DF$Effect_Ex_4
# GSE94767_se <- sqrt(GSE94767_DF$EffectVar_Ex_4)
# #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
# 
# CIT_DF <- DF[rownames(DF) == 'LAMTOR2', ]
# CIT_ES <- CIT_DF$Effect_Ex_5
# CIT_se <- sqrt(CIT_DF$EffectVar_Ex_5)
# #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
# 
# 
# ## Variable set-up
# 
# TE <- c(TCGA_ES, GSE21034_ES , GSE70768_ES, GSE94767_ES, CIT_ES )
# seTE <- c(TCGA_se,GSE21034_se, GSE70768_se, GSE94767_se, CIT_se)
# 
# # Join the variables to create a data frame
# df2_LAMTOR2 <- data.frame(GEO.study,TE,seTE )
# df2_LAMTOR2
# rownames(df2_LAMTOR2) <- df2_LAMTOR2$GEO.study
# 
# data.frame(df2_LAMTOR2, stringsAsFactors = TRUE)
# class(df2_LAMTOR2)
# str(df2_LAMTOR2)
# View(df2_LAMTOR2)
# 
# 
# m.gen <- metagen(TE = TE,
#                  seTE = seTE,
#                  studlab = GEO.study,
#                  data = df2_LAMTOR2,
#                  sm = "SMD",
#                  #method.ci = 't',
#                  method.smd = "Hedges",
#                  fixed = FALSE,
#                  random = TRUE,
#                  method.tau = "DL",#!!!!!!!!!!! here we use DJ (DerSimonian and Laird estimator) because this is what GeneMeta uses!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#                  #method.tau = "REML",
#                  hakn = FALSE, #if use false then metafor/GeneMeta/meta all matched!!!!!!!!! AND it uses Wald-type test
#                  #hakn = TRUE,
#                  title = "LAMTOR2")
# m.gen
# summary(m.gen)
# m.gen$pval.random
# 
# pdf(file = "meta_forest_LAMTOR2.pdf", width = 14, height = 8)
# forest(m.gen, studlab = TRUE)
# dev.off()





# ##################################################################################################################################################################################
# ################################################################## For DANCR #####################################################################################################
# ##################################################################################################################################################################################
# 
# # We need to create ExpressionSet
# 
# # TCGA
# eset_TCGA <- ExpressionSet(assayData = assayDataNew(exprs = as.matrix(data_TCGA)))
# pData(eset_TCGA) <-pheno_TCGA
# View(pData(eset_TCGA))
# eset_TCGA
# # eset_TCGA <- eset_TCGA[-(rownames(eset_TCGA)=="DANCR") , ]
# condition_TCGA <- as.factor(pheno_TCGA$Group)
# levels(condition_TCGA) <- c(1,0) #c(0,1)
# condition_TCGA<- as.numeric(as.character(condition_TCGA))
# eset_TCGA
# 
# # GSE21034
# 
# eset_GSE21034 <- ExpressionSet(assayData = assayDataNew(exprs = as.matrix(data_GSE21034)))
# pData(eset_GSE21034) <-pheno_GSE21034
# View(pData(eset_GSE21034))
# eset_GSE21034
# # eset_GSE21034 <- eset_GSE21034[-(rownames(eset_GSE21034)=="DANCR") , ]
# condition_GSE21034 <- as.factor(pheno_GSE21034$Group)
# levels(condition_GSE21034) <- c(1,0) #c(0,1)
# condition_GSE21034<- as.numeric(as.character(condition_GSE21034))
# eset_GSE21034
# 
# # GSE70768
# 
# eset_GSE70768 <- ExpressionSet(assayData = assayDataNew(exprs = as.matrix(data_GSE70768)))
# pData(eset_GSE70768) <-pheno_GSE70768
# View(pData(eset_GSE70768))
# eset_GSE70768
# # eset_GSE70768 <- eset_GSE70768[-(rownames(eset_GSE70768)=="DANCR") , ]
# condition_GSE70768 <- as.factor(pheno_GSE70768$Group)
# levels(condition_GSE70768) <- c(1,0) #c(0,1)
# condition_GSE70768<- as.numeric(as.character(condition_GSE70768))
# eset_GSE70768
# 
# 
# 
# 
# # # GSE94767
# # 
# # eset_GSE94767 <- ExpressionSet(assayData = assayDataNew(exprs = as.matrix(data_GSE94767)))
# # pData(eset_GSE94767) <-pheno_GSE94767
# # View(pData(eset_GSE94767))
# # eset_GSE94767
# # #eset_GSE94767 <- eset_GSE94767[-(rownames(eset_GSE94767)=="DANCR") , ] # dont need as it does have DANCR
# # condition_GSE94767 <- as.factor(pheno_GSE94767$Group)
# # levels(condition_GSE94767) <- c(1,0) #c(0,1)
# # condition_GSE94767<- as.numeric(as.character(condition_GSE94767))
# # eset_GSE94767
# 
# 
# 
# # CIT
# 
# eset_CIT <- ExpressionSet(assayData = assayDataNew(exprs = as.matrix(data_CIT)))
# pData(eset_CIT) <-pheno_CIT
# View(pData(eset_CIT))
# eset_CIT
# # eset_CIT <- eset_CIT[-(rownames(eset_CIT)=="DANCR") , ]
# condition_CIT <- as.factor(pheno_CIT$Group)
# levels(condition_CIT) <- c(1,0) #c(0,1)
# condition_CIT<- as.numeric(as.character(condition_CIT))
# eset_CIT
# 
# 
# 
# 
# 
# table(featureNames(eset_TCGA)==featureNames(eset_CIT))
# table(featureNames(eset_TCGA)==featureNames(eset_CIT))
# 
# 
# 
# ###################################################
# ### code chunk number 13: claculationAllinOne
# ###################################################
# esets     <- list(eset_TCGA,
#                   eset_GSE21034,
#                   eset_GSE70768,
#                   #eset_GSE94767,
#                   eset_CIT)
# # data.ER   <-pData(Nevins)[,"ER.status"]
# # levels(data.ER) <- c(0,1)
# # data.ER<- as.numeric(as.character(data.ER))
# classes   <- list(condition_TCGA,
#                   condition_GSE21034,
#                   condition_GSE70768,
#                   #condition_GSE94767,
#                   condition_CIT)
# #theScores <- zScores(esets,classes,useREM=FALSE,CombineExp=1:2)
# theScores <- zScores(esets,classes,useREM=TRUE,CombineExp=1:4) # CombineExp can be set to combine wanted datasets
# 
# 
# 
# 
# ###################################################
# ### code chunk number 14: show results
# ###################################################
# View(theScores)
# theScores[1:2,]
# 
# # output data
# #write.table(theScores,"Table1_Meta_synovial.txt",row.names = T,col.names = NA,quote = F,sep="\t")
# write.csv(theScores, 'Table1_Meta_DANCR.csv')# alternativelt save as csv, probably better
# save(theScores, file = "Table1_Meta_DANCR.rda") # same the unfiltered normalised expression matrix
# 
# 
# 
# # FDR
# set.seed(123)
# ScoresFDR <- zScoreFDR(esets, classes, useREM=TRUE, nperm=50,CombineExp=1:4)
# View(ScoresFDR$two.sided)
# #View(ScoresFDR$pos)
# 
# ScoresFDR_two.sided <- ScoresFDR$two.sided
# View(ScoresFDR_two.sided)
# write.csv(ScoresFDR_two.sided, 'Table2_Meta_DANCR_FDR_default.csv')# alternativelt save as csv, probably better
# save(ScoresFDR_two.sided, file = "Table2_Meta_DANCR_FDR_default.rda") # same the unfiltered normalised expression matrix
# 
# 
# 
# # BH-FDR
# 
# ## https://www.youtube.com/watch?v=DRRT4k7A2dM
# 
# b<- as.data.frame(ScoresFDR_two.sided)
# class(b)
# b$rawPbyZ <- 2*pnorm(-abs(b$zSco))
# View(b)
# b$FDR_BH <- p.adjust(b$rawPbyZ, "BH") # I dont think it is necessary to do fdr if we assuming to looking gene one by one
# View(b)
# write.csv(b, 'Table3_Meta_DANCR_FDR_BH.csv')# alternativelt save as csv, probably better
# 
# 
# 
# ###################################################
# ### Pooled ES volcano
# ###################################################
# 
# GeneMeta_res<- as.data.frame(ScoresFDR_two.sided)
# GeneMeta_res$pooled_ES <- GeneMeta_res$MUvals
# View(GeneMeta_res)
# 
# deg.data <- as.data.frame(GeneMeta_res[ , c(17,16)])
# View(deg.data)
# class(deg.data)
# deg.data$logP <- -log10(deg.data$Chisq)
# 
# deg.data$Group = "Not significant"
# deg.data$Group[which( (deg.data$Chisq < 0.05))] = "Significant"
# #deg.data$Group[which( (deg.data$adj.P.Val < p.adjust.cutoff) & (deg.data$logFC < -log(fold.change.cutoff,2)) )] = "Down-regulated"
# deg.data$Label = ""
# deg.data$LabelSF = ""
# deg.data <- deg.data[order(deg.data$Chisq), ]
# up.genes <- head(rownames(deg.data)[which(deg.data$Group == "Significant")], 10)
# #down.genes <- head(deg.data$Symbol[which(deg.data$Group == "Down-regulated")], 10)
# deg.top10.genes <- (as.character(up.genes))
# deg.data$Label[match(deg.top10.genes, rownames(deg.data))] <- deg.top10.genes
# deg.data$LabelSF[match(deg.top10.genes, rownames(deg.data))] <- deg.top10.genes # NEED to run ' Extracting splicing factor info' below beforehand for 'degSF'
# deg.data$Label_all <- rownames(deg.data)
# View(deg.data)
# 
# library(ggplot2)
# library(ggpubr)
# p1 <- ggscatter(deg.data, x = "pooled_ES", y = "logP",
#                 color = "Group",
#                 shape = "Group",
#                 palette = c("#616161", "#000000"),#c("#BBBBBB", "#CC0000"),
#                 size = 2.5,
#                 label = deg.data$Label_all,
#                 #label = deg.data$LabelSF,
#                 font.label = 12,
#                 repel = T,
#                 xlab = "Pooled effect size (Hedges'g)",
#                 ylab = "-log10(p-value)",
#                 ylim = c(0,16)) +
#   theme_bw() +
#   geom_hline(yintercept = -log(0.05,10), linetype="dashed") +
#   ggtitle("Volcano plot of selected gene markers")
# #geom_vline(xintercept = c(log(fold.change.cutoff,2),-log(fold.change.cutoff,2)), linetype="dashed")+
# #ggtitle(paste("Volcano plot","(",colnames(fit2$contrasts),")"))
# p1
# library(cowplot)
# ## output plot
# save_plot("volcano_fourSets_DANCR.pdf",p1,base_width = 6, base_height = 6)
# save_plot("volcano_fourSets_DANCR.tiff",p1,base_width = 6, base_height = 6,
#           compression = "lzw",dpi = 500)
# 
# 
# 
# 
# #-----------------------------------------------------------------##
# #                   Forest for DANCR only                         ##
# #-----------------------------------------------------------------##   
# library(meta)
# 
# 
# #https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/pooling-es.html
# #https://bookdown.org/MathiasHarrer/Doing_Meta_Analysis_in_R/effects.html
# 
# lname <- load('Table1_Meta_DANCR.rda')
# lname
# 
# 
# DF <- as.data.frame(theScores)
# 
# class(DF)
# 
# 
# View(DF)
# 
# 
# ## DANCR
# GEO.study <- c('TCGA','GSE21034','GSE70768', 'CIT')
# 
# TCGA_DF <- DF[rownames(DF) == 'DANCR', ]
# TCGA_ES <- TCGA_DF$Effect_Ex_1
# TCGA_se <- sqrt(TCGA_DF$EffectVar_Ex_1)
# #GSE57218_se <- sqrt(((5+5)/(5*5))+((GSE57218_ES)^2/(2*(5+5))))
# 
# GSE21034_DF <- DF[rownames(DF) == 'DANCR', ]
# GSE21034_ES <- GSE21034_DF$Effect_Ex_2
# GSE21034_se <- sqrt(GSE21034_DF$EffectVar_Ex_2)
# #GSE117999_se <- sqrt(((10+10)/(10*10))+((GSE117999_ES)^2/(2*(10+10))))
# 
# 
# GSE70768_DF <- DF[rownames(DF) == 'DANCR', ]
# GSE70768_ES <- GSE70768_DF$Effect_Ex_3
# GSE70768_se <- sqrt(GSE70768_DF$EffectVar_Ex_3)
# #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
# 
# # GSE94767_DF <- DF[rownames(DF) == 'DANCR', ]
# # GSE94767_ES <- GSE94767_DF$Effect_Ex_4
# # GSE94767_se <- sqrt(GSE94767_DF$EffectVar_Ex_4)
# # #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
# 
# CIT_DF <- DF[rownames(DF) == 'DANCR', ]
# CIT_ES <- CIT_DF$Effect_Ex_4
# CIT_se <- sqrt(CIT_DF$EffectVar_Ex_4)
# #GSE169077_se <- sqrt(((10+9)/(10*9))+((GSE169077_ES)^2/(2*(10+9))))
# 
# 
# ## Variable set-up
# 
# TE <- c(TCGA_ES, GSE21034_ES , GSE70768_ES, CIT_ES )
# seTE <- c(TCGA_se,GSE21034_se, GSE70768_se, CIT_se)
# 
# # Join the variables to create a data frame
# df2_DANCR <- data.frame(GEO.study,TE,seTE )
# df2_DANCR
# rownames(df2_DANCR) <- df2_DANCR$GEO.study
# 
# data.frame(df2_DANCR, stringsAsFactors = TRUE)
# class(df2_DANCR)
# str(df2_DANCR)
# View(df2_DANCR)
# 
# 
# m.gen <- metagen(TE = TE,
#                  seTE = seTE,
#                  studlab = GEO.study,
#                  data = df2_DANCR,
#                  sm = "SMD",
#                  #method.ci = 't',
#                  method.smd = "Hedges",
#                  fixed = FALSE,
#                  random = TRUE,
#                  method.tau = "DL",#!!!!!!!!!!! here we use DJ (DerSimonian and Laird estimator) because this is what GeneMeta uses!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#                  #method.tau = "REML",
#                  hakn = FALSE, #if use false then metafor/GeneMeta/meta all matched!!!!!!!!! AND it uses Wald-type test
#                  #hakn = TRUE,
#                  title = "DANCR")
# m.gen
# summary(m.gen)
# m.gen$pval.random
# 
# pdf(file = "meta_forest_DANCR.pdf", width = 14, height = 8)
# forest(m.gen, studlab = TRUE)
# dev.off()







