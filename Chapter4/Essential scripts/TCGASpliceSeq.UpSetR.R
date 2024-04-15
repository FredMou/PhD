



#======================================= UpSet plot for all TCGA SpliceSeq events
## Refer to 'DEAS_final' in '129immAS -> workspace' folder
library(limma)
library(impute) 
library(UpSetR)           
asFile="asMatrix.txt"     
setwd("/Users/zhuofanmou/Documents/PhD Exeter/ClariomD/alternative_splicing/R_workspace_v2/UpsetR")  # Working directory


rt=read.table("asMatrix.txt", header=T, sep="\t", check.names=F, row.names=1)
View(rt)

# UpSet
gene=sapply(strsplit(rownames(rt),"\\|"), "[", 1)
asType=sapply(strsplit(rownames(rt),"\\|"), "[", 3)
upsetList=list(AA=unique(gene[asType=="AA"]),
               AD=unique(gene[asType=="AD"]),
               AP=unique(gene[asType=="AP"]),
               AT=unique(gene[asType=="AT"]),
               ES=unique(gene[asType=="ES"]),
               ME=unique(gene[asType=="ME"]),
               RI=unique(gene[asType=="RI"]) )
upsetData=fromList(upsetList)

#plot and save
pdf(file="upset_allEvents_TCGASpliseseq.pdf", width=12, height=6, onefile=FALSE)
upset(upsetData,
      nsets = 7,                   # number of AS event type
      nintersects = 50,            # number of gene intersection set blogs/bars
      order.by = "freq",           # sort by number of intersected genes in each blog/bar
      show.numbers = "yes",        
      number.angles = 20,          
      point.size = 1.5,           
      matrix.color="red",         
      line.size = 0.8,             
      mainbar.y.label = "Gene Intersections", 
      sets.x.label = "Set Size")
dev.off()