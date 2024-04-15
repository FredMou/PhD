# working directory set up
setwd('/Users/zhuofanmou/Documents/PhD Exeter/ClariomD/alternative_splicing/R_workspace_v2/limmaVoom_DGEA_TCGA')
getwd()






options(stringsAsFactors = F)
library(stringr)
library(limma)
library(cowplot)
library(openxlsx)
library(ggpubr)
library(ggthemes)

# read in expression value
data = read.table("datExp.normalized_filtered.txt",header = T,sep="\t",row.names = 1,check.names = F)
View((data))

# read in phenotype
pheno0 = read.table("target.txt",header = T,sep="\t")
View((pheno0))
#pheno0[,2] = str_remove_all(pheno0[,2]," +$")

pheno = pheno0[match(colnames(data),pheno0[,1]),]
table(colnames(data) == pheno$barcode)
View((pheno))
#pheno$Group[pheno$Group == "Healthy"] <- "normal"

# 
# 
# # assess for potential covariates----
# str(pheno)
# pheno$Sex <- factor(pheno$Sex)
# pheno$Group2 <- factor(pheno$Group)
# pheno$Group2 <- as.numeric(pheno$Group2)
# pheno$Group2[pheno$Group2 == 1] <- 0
# pheno$Group2[pheno$Group2 == 2] <- 1
# #pheno$Group2 <- as.numeric(pheno$Group2)  # as factor or as numerics won't change the results!!!
# pheno$Group2 <- as.factor(pheno$Group2)
# View((pheno))
# str(pheno)
# 
# 
# ## for age (significant)
# age <- glm(pheno$Group2 ~ pheno$Age, family = binomial)
# summary(age)
# ## for sex (not significant)
# sex <- chisq.test(table(pheno$Group2, pheno$Sex))
# sex <- fisher.test(table(pheno$Group2, pheno$Sex))
# (sex)
# 
# ## limma
# 
# # potential covariates
# 
# str(pheno)
#pheno$Sex <- factor(pheno$Sex)
#model design
Group = factor(pheno$Group,levels=c('Tumour','Normal'))
design = model.matrix(~0+Group)
colnames(design) <- c('Tumour','Normal')

# design = model.matrix(~0+Group+pheno$Sex+pheno$Age)
# colnames(design) <- c('OA','normal','Sex','Age')
design


#linear model fitness
fit <- lmFit(data, design)

#generate contrast matrix
contrast.matrix <- makeContrasts(Tumour-Normal, #1; logFC positive=up-regulated between the group£»logFC negative=down-regulated between the group
                                 levels=design)
#constrast model fit 
fit2 <- contrasts.fit(fit, contrast.matrix)

#bayes model 
fit2 <- eBayes(fit2)

#get DEGs
p.adjust.cutoff = 0.05
fold.change.cutoff = 2 #no change; other (stringent) options:4/2/1.5

diff = topTable(fit2,adjust.method="fdr",
                p.value=p.adjust.cutoff,
                lfc=log(fold.change.cutoff,2),
                number=99999,sort.by = 'logFC')
View(diff)

#add annotations
#anno = read.csv("Annotation.csv",header = T, na.strings = c("","NA")) # Gene symbol is not label as 'NA' but "", so we need 'na.strings = c("","NA")'
anno = read.csv("Annotation.csv",header = T) # Gene symbol is not label as 'NA' but "", so we need 'na.strings = c("","NA")'

anno<- sapply(anno, as.character)
anno[is.na(anno)] <- ""
anno <- as.data.frame(anno)
class(anno)
str(anno)
View((anno))

colnames(anno) = c("ID","Gene.Symbol")
View((anno))

diff$ID = rownames(diff)
diff$Gene = anno$Gene.Symbol[match(diff$ID,anno$ID)] #  (Based on the 'diff' results) we want to MATCH/FIND the 'position/indenx' of the 'diff' rownames/probIDs IN the anno$ID (presented as in probeID). because anno$Gene.Symbol and anno$ID are in the same order, we then use the index to find the correct position of the rows for gene symbol
#diff$EntrezID = anno$EntrezID[match(diff$ID,anno$ID)] #  (Based on the 'diff' results) we want to MATCH/FIND the 'position/indenx' of the 'diff' rownames/probIDs IN the anno$ID (presented as in probeID). because anno$Gene.Symbol and anno$ID are in the same order, we then use the index to find the correct position of the rows for gene symbol
#diff$EntrezID <- as.character(diff$EntrezID)
diff$Compare = colnames(fit2$contrasts)
str(diff)
View(diff)
table(duplicated(diff$Gene))
# filteration
diff = diff[,c(7,8,1:6,9)] # if AveExpr column appears in the 'diff', depends on the 'limma' version; 3.46.0 has the column BUT NOT in 3.50.1 
diff = diff[diff$Gene != "---",] # remove probes without gene symbols (some GEO use '---')
# diff$Gene <- gsub("MIR636 /// SRSF2", "SRSF2", diff$Gene)  # replace one-to-many by selecting the SF one (probably not a good idea......);  SPACE!!!!!!!!
#diff = diff[!str_detect(diff$Gene,"///"),] # remove the rest probes matching to multiple genes
diff = diff[diff$Gene != "" , ] #!!!remove probes without gene symbols (some GEO use blank '' NOT even 'NA')
View((diff))
dim(diff) #3696*10
# output DEGs
write.xlsx(diff,"Table1.DEG.xlsx")




##----------------------------------------##
##  specification for interested gene set ##
##----------------------------------------##

#genelist <- data.frame(id=c("MKI67","PCNA","BCL2","GPX4","TFRC"),genelist=c("MKI67","PCNA","BCL2","GPX4","TFRC"))
SFlist = read.table("SFlist.txt",header = F,sep="\t")
class(SFlist)
View(SFlist)
colnames(SFlist) <- SFlist[1,]
SFlist <- SFlist[-1,]


# all genes with limma results

## 1 For meta-analysis
limma_meta = topTable(fit2,adjust.method="fdr",coef=1,p.value=1,
                      lfc=log(1,2),confint = TRUE,number=99999,sort.by = 'p') #  include all genes; here we DO NOT want to sort/rank the list by any column variables because later we want to match and extract the indecies of the genes (by gene symbol from deg.data) that we are interested between the probe ID of 'fit2' and 'deg.data'. if we sort the order here, then the order of probe ID between the 'fit2' and 'deg.data' will be DIFFERENT!!!

limma_meta$ID = rownames(limma_meta)
limma_meta$Gene = anno$Gene.Symbol[match(limma_meta$ID , anno$ID)]
limma_meta$EntrezID = anno$EntrezID[match(limma_meta$ID , anno$ID)]
limma_meta$Compare = colnames(fit2$contrasts)
limma_meta = limma_meta[ , c(9,10,1:8,11)]
limma_meta = limma_meta[limma_meta$Gene != "---",]
#limma_meta$Gene <- gsub("MIR636 /// SRSF2", "SRSF2", limma_meta$Gene)
#limma_meta = limma_meta[!str_detect(limma_meta$Gene,"///"),]
limma_meta = limma_meta[limma_meta$Gene != "" , ]
dim(limma_meta) #13227*12
View(limma_meta)

#!!! when the are multiple probes to the same gene, we take the gene/probe that has the highest average across the sample!!!
limma_meta <- limma_meta %>%
  arrange(Gene , -AveExpr) %>%
  filter(!duplicated(Gene))
dim(limma_meta) # 13339*12
View(limma_meta)


# Check if there are still duplicates
table(duplicated(limma_meta$Gene))
# output all probes from limma
write.xlsx(limma_meta,"Table3.allGenes.limma.meta.xlsx")




## 2 for this script and following work
deg.data = topTable(fit2,adjust.method="fdr",coef=1,p.value=1,
                    lfc=log(1,2),number=99999,sort.by = 'none') #  include all genes; here we DO NOT want to sort/rank the list by any column variables because later we want to match and extract the indecies of the genes (by gene symbol from deg.data) that we are interested between the probe ID of 'fit2' and 'deg.data'. if we sort the order here, then the order of probe ID between the 'fit2' and 'deg.data' will be DIFFERENT!!!

View(deg.data)

# deg.data = topTable(fit2,adjust.method="fdr",coef=1,p.value=1,
#                     lfc=log(fold.change.cutoff,2),number=50000,sort.by = 'logFC')


# Match the corresponding gene symbols
deg.data$Symbol = anno$Gene.Symbol[match(rownames(deg.data),anno$ID)]
deg.data$EnsembID = anno$ID[match(rownames(deg.data),anno$ID)]

# #check if there is a one-to-multiple case
# check <- deg.data[deg.data$Symbol == "MIR636 /// SRSF2",]
# check

# #pick one gene of interest if one-to-multiple case occur (probably not a good idea)
# deg.data$Symbol <- gsub("MIR636 /// SRSF2", "SRSF2", deg.data$Symbol)
# check2 <- deg.data[deg.data$Symbol == "SRSF2",]
# check2


# Extract the limma results for the interested SF genes
SF_limma <- subset(deg.data, EnsembID %in% SFlist$ENSEMBL)
View(SF_limma)


# diff_SFlist1 <- deg.data[match(SFlist$V1, deg.data$Symbol),]
# View(diff_SFlist1)

# only apply multiple tests on the interested SF list!!!----

# If want to obtain all probesID (regardeless duplicated gene symbols); matching by the probeID!!!

allSF_fdr_probeID = topTable(fit2[ match(rownames(SF_limma), rownames(deg.data)), ],adjust.method="fdr",
                             p.value=1,
                             lfc=log(1,2),
                             number=99999,sort.by = "logFC") # sort.by here WILL NOT affect the adj.p calculation like 'BH',  as it is done after the multiple tests; so depending on the adjust.method we pick, it will do their job automatically (e.g., by ranking the raw p values using 'BH') and DO NOT need to manually ranking it beforehand
View(allSF_fdr_probeID)


# combine the gene symbol in SF_limma and all_SF_probeID according to the probeID

probeID_interested <- as.data.frame(SF_limma[ , -c(1:5) ])
probeID_interested$Gene_symbol <- probeID_interested$Symbol
probeID_interested <- as.data.frame(probeID_interested[ , -1 ])
View(probeID_interested)


same=intersect(row.names(probeID_interested), row.names(allSF_fdr_probeID))
length(probeID_interested)
df1=probeID_interested[same,]
head(rownames(df1))
df2=allSF_fdr_probeID[same,]
head(rownames(df2))
SF_probes_final=cbind(df1,df2)
SF_probes_final <- SF_probes_final[ , -1]
SF_probes_final$ProbeID <- rownames(SF_probes_final)
SF_probes_final <- SF_probes_final[ , c(9,1:8)]
View(SF_probes_final)





# # If want to obtain under certain cut-offs
# cutSF_fdr_probeID = topTable(fit2[ match(rownames(SF_limma), rownames(deg.data)), ],adjust.method="fdr",
#                         p.value=p.adjust.cutoff,
#                         lfc=log(1,2),
#                         number=50000,sort.by = "logFC") # sort.by here WILL NOT affect the adj.p calculation like 'BH',  as it is done after the multiple tests; so depending on the adjust.method we pick, it will do their job automatically (e.g., by ranking the raw p values using 'BH') and DO NOT need to manually ranking it beforehand
# 
# 
# View(cutSF_fdr_probeID) # adj.p from it are not the same as in the 'diff_SFlist1', which computed from the full data

# (checked) ALTERNATIVE way of calculating the BH adj p values, easier when we have a specific genes of interest and want to calculate adj p for that set ONLY
# NOTE! This is ONLY applicable when for the total genes of interest (i.e., unfiltered under the P&log2FC cut-offs above), otherwise, the adj.p will be different
check3 <- SF_probes_final
check3$fdr <- p.adjust(SF_probes_final$P.Value, "BH")
View(check3)

# output SF DEGs
write.xlsx(SF_probes_final,"Table2.SF.DEG.xlsx")




##----------------------------------------##
##    hierachical clustering heatmap2      ##
##----------------------------------------##
###hierachical clustering heatmap for significant DE probes (splicing factors) only----
SF_probes_final_fdr0.05 <- subset(SF_probes_final , adj.P.Val <= 0.05)
View(SF_probes_final_fdr0.05) 
dim(SF_probes_final_fdr0.05) # 49 DESFs
table(duplicated(SF_probes_final_fdr0.05$EnsembID))
table(duplicated(SF_probes_final_fdr0.05$Gene_symbol))

# get expression values of significant DEG SF probes
degSFprobes.exp = data[match(rownames(SF_probes_final_fdr0.05),rownames(data)),]
degSFprobes.exp$probeID = anno$ID[match(rownames(degSFprobes.exp),anno$ID)]

degSFprobes.exp$geneSymbol = SF_probes_final_fdr0.05$Gene_symbol[match(rownames(degSFprobes.exp),SF_probes_final_fdr0.05$ProbeID)]

degSFprobes.exp = degSFprobes.exp[!duplicated(degSFprobes.exp$probeID),]
rownames(degSFprobes.exp) = degSFprobes.exp$probeID
library(stringr)
degSFprobes.exp$ID <- str_c(degSFprobes.exp$probeID , ' | ' ,degSFprobes.exp$geneSymbol)

#rownames(degSFprobes.exp) = degSFprobes.exp$ID # NO NEED for TCGA fata
rownames(degSFprobes.exp) = degSFprobes.exp$geneSymbol

degSFprobes.exp = degSFprobes.exp[,1:(ncol(degSFprobes.exp) - 3)]
View(degSFprobes.exp)

deSFexp_forOmicShare <- rbind(ID=colnames(degSFprobes.exp),degSFprobes.exp)
View(deSFexp_forOmicShare)
write.table(deSFexp_forOmicShare,file="Table4.sigSF.txt",sep="\t",col.names=F,quote=F)

# output SF Expression
#write.xlsx(deSFexp_forOmicShare,"Table4.sigSF.exp.xlsx")

## draw heatmap
library(pheatmap)
annotation_col = as.data.frame(pheno$Group)
colnames(annotation_col) = "Group"
#rownames(annotation_col) = colnames(degSFprobes.exp)
rownames(annotation_col) = pheno$sample
color.key <- c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000")
sample.color.list = list(Group=c(get_palette("nejm",length(unique(annotation_col[,1])))))
names(sample.color.list$Group) = unique(annotation_col[,1])
p2 = pheatmap(degSFprobes.exp, color = colorRampPalette(color.key)(50), 
              border_color = NA,
              annotation_col = annotation_col,
              labels_row = NULL,clustering_method = "ward.D2",
              show_rownames = T,show_colnames = T,fontsize_col = 10,
              annotation_colors = sample.color.list
)



## output plot 
save_plot("Figure2B.heatmap.pdf",p2,base_width = 8, base_height = 6)
save_plot("Figure2B.heatmap.tiff",p2,base_width = 8, base_height = 6,
          compression = "lzw",dpi = 500)





