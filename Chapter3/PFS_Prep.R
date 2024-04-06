#install.packages("survival")
#install.packages("survminer")


#Setup----
library(survival)
library(survminer)

cliFile="Survival_SupplementalTable_S1_20171025_xena_sp"     #Clinical data import (including PFS)
setwd("/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Model construction/R_data_needed")     #???ù???Ŀ¼

#??ȡ?ٴ??????ļ?----
cli_PAN=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
View(cli_PAN)
dim(cli_PAN)
# cli_FH <- read.table(cliFile2, header=T, sep="\t", check.names=F, row.names=1)
# dim(cli_FH)
# View(cli_FH)

#From UCSC PANCANCER to obtain PRAD clinical info----
cli_PAN_PFS=cli_PAN[,c("PFI.time","PFI")]
cli_PAN_DFS=cli_PAN[,c("DFI.time","DFI")]
#cli_PAN_DSS=cli_PAN[,c("DSS.time","DSS")]
cli_PAN_PFS=na.omit(cli_PAN_PFS)
cli_PAN_DFS=na.omit(cli_PAN_DFS)
#cli_PAN_DFS=na.omit(cli_PAN_DSS)
dim(cli_PAN_PFS)
dim(cli_PAN_DFS)
#dim(cli_PAN_DSS)
View(cli_PAN_PFS)
View(cli_PAN_DFS)
#View(cli_PAN_DSS)
colnames(cli_PAN_PFS)=c('days_to_last_follow_up','vital_status')
colnames(cli_PAN_DFS)=c('days_to_last_follow_up','vital_status')
#colnames(cli_PAN_DSS)=c('days_to_last_follow_up','vital_status')
any(duplicated(rownames(cli_PAN_PFS))) # check if there is duplicated names
any(duplicated(rownames(cli_PAN_DFS))) # check if there is duplicated names
#any(duplicated(rownames(cli_PAN_DSS))) # check if there is duplicated names


##EdgeR
lnames_edgeR <- load('/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Model construction/R_data_needed/ExpDat_norm_log2_edgeR.rda')
lnames_edgeR

View(datNorm_all_log2_edgeR)
exp <- as.data.frame(datNorm_all_log2_edgeR)
exp <- exp[ , c('HNRNPL', 'ELAVL1', 'PCBP1', 'PCBP2', 'PABPN1','PTPRF','MRPL24','DANCR', 'MYC', 'TRPM4')]
View(exp)
dim(exp)



##??????????1:????????A??Bȥ???? ??????PANCANCER; edgeR----
d<-exp
d$names <- rownames(d)
View(d)
library(stringr)
library(dplyr)
d$names = gsub(pattern = "-11R.*", replacement = "", x = d$names)
d$names = gsub(pattern = "-22R.*", replacement = "", x = d$names)
d$names = gsub(pattern = "-21R.*", replacement = "", x = d$names)
d$names = gsub(pattern = "-12R.*", replacement = "", x = d$names)
d$names = gsub(pattern = "-31R.*", replacement = "", x = d$names)
d$names = gsub(pattern = "-04R.*", replacement = "", x = d$names)
d$names = gsub(pattern = "-41R.*", replacement = "", x = d$names)
d$names = gsub(pattern = "-05R.*", replacement = "", x = d$names)
d$names = gsub(pattern = "-01R.*", replacement = "", x = d$names)
d$names = gsub(pattern = "-13R.*", replacement = "", x = d$names)
d$names = gsub(pattern = "-32R.*", replacement = "", x = d$names)
d$names = gsub(pattern = "-61R.*", replacement = "", x = d$names)
d$names = gsub(pattern = "-02R.*", replacement = "", x = d$names)

d <- d %>%
  mutate_at("names", str_replace, "-11A", "-11")
d <- d %>%
  mutate_at("names", str_replace, "-11B", "-11")
d <- d %>%
  mutate_at("names", str_replace, "-01A", "-01")
d <- d %>%
  mutate_at("names", str_replace, "-01B", "-01")

View(as.matrix(rownames(d)))
View(as.matrix(d$names))
d <- d[!duplicated(d$names),]
rownames(d) <- d$names
any(duplicated(d$names)) # check if there is duplicated names
table(duplicated(d$names)) 
View(d)
class(d)
d <- d[!grepl("-11", d$names),]
dim(d)
d <- d[ , -ncol(d)]
exp_edgeR_TCGA <- d
View(exp_edgeR_TCGA)
dim(exp_edgeR_TCGA)







#?ϲ?????, PANCANCER PFS edgeR----
sameSample_PAN_PFS_edgeR=intersect(row.names(cli_PAN_PFS), row.names(exp_edgeR_TCGA))
length(sameSample_PAN_PFS_edgeR)
cli_PAN_PFS_edgeR=cli_PAN_PFS[sameSample_PAN_PFS_edgeR,]
head(rownames(cli_PAN_PFS_edgeR))
exp_PFS_edgeR=exp_edgeR_TCGA[sameSample_PAN_PFS_edgeR,]
head(rownames(exp_PFS_edgeR))
rt_final_PFS_edgeR=cbind(cli_PAN_PFS_edgeR, exp_PFS_edgeR)
View(rt_final_PFS_edgeR)
rt_final_PFS_edgeR$days_to_last_follow_up=(rt_final_PFS_edgeR$days_to_last_follow_up)/365
table(rt_final_PFS_edgeR$vital_status)
dim(rt_final_PFS_edgeR)
any(duplicated(rownames(rt_final_PFS_edgeR))) # check if there is duplicated names
save(rt_final_PFS_edgeR, file = "TCGA_PAN_PFS_GSig_Input_edgeR.rda") # Final PFS+10 gene expression dataset











