# setwd("/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/Boxplot_expression")
# getwd()

#Refer to KTM '2.AdvancedVis.R' immune infrastration

##EdgeR expression
lnames_edgeR <- load('/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Model construction/R_data_needed/ExpDat_norm_log2_edgeR.rda')
lnames_edgeR

View(datNorm_all_log2_edgeR)
exp <- as.data.frame(datNorm_all_log2_edgeR)
exp <- exp[ , c('HNRNPL', 'ELAVL1', 'PCBP1', 'PCBP2', 'PABPN1','PTPRF', 'MRPL24','DANCR', 'MYC', 'TRPM4')]
View(exp)
dim(exp)


# Group for Normal vs tumour
data <- t(exp)
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data_modi=rbind(data, Type=as.numeric(group))
data_modi=as.data.frame(t(data_modi[c(rownames(data_modi),"Type"),]))
data_modi$Type=ifelse(data_modi$Type==0,"Tumor","Normal")
dim(data_modi)
head(data_modi)
data_modi <- data_modi[ , -ncol(data_modi)]
View(data_modi)

geneID <- colnames(exp)
geneID

# Normal vs tumour boxplot
plot.info <- NULL
for (i in 1:length(geneID)) {
  idx.sub <- which(colnames(data_modi) == geneID[i])
  sub <- data.frame(
    PATIENT_ID = rownames(data_modi),
    Gene_symbol = geneID[i],
    Group = data_modi$Type,
    Expression = data_modi[, idx.sub]
  )
  plot.info <- rbind(plot.info, sub)
}

library(ggplot2)
library(tidyverse)
library(ggpubr)
# box plot group  notmal vs tumour
pdf(file="Box_TvN_10.pdf",height=6,width=10)

p1_statusComp <- ggboxplot(plot.info, 
                x = "Gene_symbol", y = "Expression",
                color = "black", fill = "Group",palette = "lancet",
                xlab = "", ylab = "Gene expression",
                main = "Gene expression comparison (normal vs. tumour)")+ #between normal (52) and tumour (499)
  stat_compare_means(aes(group = Group), label = "p.signif", method = "wilcox.test", hide.ns = F)+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30,hjust = 1))
  
  # stat_compare_means(label = "p.signif", method = "wilcox.test", # wilcox.test
  #                    ref.group = ".all.", hide.ns = F) +
  # theme_bw()+
  # theme(panel.background = element_blank(),
  #       panel.grid = element_blank(),
  #       axis.text.x = element_text(angle = 60,hjust = 1))
p1_statusComp

dev.off()


clinical <- read.table("/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Model construction/R_data_needed/df_clinic.genemarkers.COX_final.txt", header = T, sep = "\t") 
View(clinical)
rownames(clinical) <- clinical$barcode


#Merged data----
sameSample_expClinic=intersect(row.names(data_modi), row.names(clinical))
length(sameSample_expClinic)
clinical2=clinical[sameSample_expClinic,]
head(rownames(clinical2))
data_modi2=data_modi[sameSample_expClinic,]
head(rownames(data_modi2))
data_combined=cbind(clinical2, data_modi2)
View(data_combined)


# Gleason boxplot:  NOTE!!! dont use '≥' sign, warning will come up no matter which format of the variable (factor, numeric or character) used

library(tidyr)
data_combined_Gleason <- data_combined %>% drop_na(gleason_score)
dim(data_combined_Gleason)
View(data_combined_Gleason)
str(data_combined_Gleason)



data_combined_Gleason[data_combined_Gleason$gleason_score == '<8',]$gleason_score <- "<8"#"less than 8" # OR:0
data_combined_Gleason[data_combined_Gleason$gleason_score == '≥8',]$gleason_score <- ">=8"#"greater or equal to 8" # OR:1

data_combined_Gleason$gleason_score <- as.factor(data_combined_Gleason$gleason_score)



plot.info.Gleason <- NULL
for (i in 1:length(geneID)) {
  idx.sub.Gleason <- which(colnames(data_combined_Gleason) == geneID[i])
  sub.Gleason <- data.frame(
    PATIENT_ID = rownames(data_combined_Gleason),
    Gene_symbol = geneID[i],
    Group = data_combined_Gleason$gleason_score,
    Expression = data_combined_Gleason[, idx.sub.Gleason]
  )
  plot.info.Gleason <- rbind(plot.info.Gleason, sub.Gleason)
}


# box plot Gleason
pdf(file="Box_Gleason_10.pdf",height=6,width=10)
p2_gleasonComp <- ggboxplot(plot.info.Gleason, 
                           x = "Gene_symbol", y = "Expression",
                           color = "black", fill = "Group",palette = "lancet",
                           xlab = "", ylab = "Gene expression",
                           main = "Gene expression comparison (Gleason score: <8 vs. >=8)")+  #Gleason score <8 (293) and ??8 (206)
  
  stat_compare_means(aes(group = Group), label = "p.signif", method = "wilcox.test", hide.ns = F)+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30,hjust = 1))
  # stat_compare_means(label = "p.signif", method = "wilcox.test", # wilcox.test
  #                    ref.group = ".all.", hide.ns = F) +
  # theme_bw()+
  # theme(panel.background = element_blank(),
  #       panel.grid = element_blank(),
  #       axis.text.x = element_text(angle = 60,hjust = 1))
p2_gleasonComp
dev.off()


# T stage boxplot
library(tidyr)
data_combined_T <- data_combined %>% drop_na(ajcc_pathologic_t)
dim(data_combined_T)
View(data_combined_T)

plot.info.T <- NULL
for (i in 1:length(geneID)) {
  idx.sub.T <- which(colnames(data_combined_T) == geneID[i])
  sub.T <- data.frame(
    PATIENT_ID = rownames(data_combined_T),
    Gene_symbol = geneID[i],
    Group = data_combined_T$ajcc_pathologic_t,
    Expression = data_combined_T[, idx.sub.T]
  )
  plot.info.T <- rbind(plot.info.T, sub.T)
}


# box plot T

pdf(file="Box_Tstage_10.pdf",height=6,width=10)
p3_TComp <- ggboxplot(plot.info.T, 
                           x = "Gene_symbol", y = "Expression",
                           color = "black", fill = "Group",palette = "lancet",
                           xlab = "", ylab = "Gene expression",
                           main = "Gene expression comparison (tumour stage: T2 vs. T3-T4)")+ # tumour stage T2 (189)  and T3-T4 (303)
  stat_compare_means(aes(group = Group), label = "p.signif", method = "wilcox.test", hide.ns = F)+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30,hjust = 1))
  
   # stat_compare_means(label = "p.signif", method = "wilcox.test", # wilcox.test
  #                    ref.group = ".all.", hide.ns = F) +
  # theme_bw()+
  # theme(panel.background = element_blank(),
  #       panel.grid = element_blank(),
  #       axis.text.x = element_text(angle = 60,hjust = 1))
p3_TComp
dev.off()

# N stage boxplot
library(tidyr)
data_combined_N <- data_combined %>% drop_na(ajcc_pathologic_n)
dim(data_combined_N)
View(data_combined_N)

plot.info.N <- NULL
for (i in 1:length(geneID)) {
  idx.sub.N <- which(colnames(data_combined_N) == geneID[i])
  sub.N <- data.frame(
    PATIENT_ID = rownames(data_combined_N),
    Gene_symbol = geneID[i],
    Group = data_combined_N$ajcc_pathologic_n,
    Expression = data_combined_N[, idx.sub.N]
  )
  plot.info.N <- rbind(plot.info.N, sub.N)
}


# box plot N
pdf(file="Box_Nstage_10.pdf",height=6,width=10)
p4_NComp <- ggboxplot(plot.info.N, 
                           x = "Gene_symbol", y = "Expression",
                           color = "black", fill = "Group",palette = "lancet",
                           xlab = "", ylab = "Gene expression",
                           main = "Gene expression comparison (lymph node status: N0 vs. N1)")+ # Lymph node status: N0(347) vs. N1(79)
  stat_compare_means(aes(group = Group), label = "p.signif", method = "wilcox.test", hide.ns = F)+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30,hjust = 1))
  
  
  # stat_compare_means(label = "p.signif", method = "wilcox.test", # wilcox.test
  #                    ref.group = ".all.", hide.ns = F) +
  # theme_bw()+
  # theme(panel.background = element_blank(),
  #       panel.grid = element_blank(),
  #       axis.text.x = element_text(angle = 60,hjust = 1))
p4_NComp
dev.off()



# Age boxplot: NOTE!!! dont use '≥' sign, warning will come up no matter which format of the variable (factor, numeric or character) used
library(tidyr)
data_combined_age <- data_combined %>% drop_na(age_at_index)
dim(data_combined_age)
View(data_combined_age)
str(data_combined_age)



data_combined_age[data_combined_age$age_at_index == '<60',]$age_at_index <- "<60"#"less than 60 years" # OR:0
data_combined_age[data_combined_age$age_at_index == '≥60',]$age_at_index <- ">=60"#"greater or equal to 60 years" # OR:1

data_combined_age$age_at_index <- as.factor(data_combined_age$age_at_index)

plot.info.age <- NULL
for (i in 1:length(geneID)) {
  idx.sub.age <- which(colnames(data_combined_age) == geneID[i])
  sub.age <- data.frame(
    PATIENT_ID = rownames(data_combined_age),
    Gene_symbol = geneID[i],
    Group = data_combined_age$age_at_index,
    Expression = data_combined_age[, idx.sub.age]
  )
  plot.info.age <- rbind(plot.info.age, sub.age)
}


# box plot Age
pdf(file="Box_age_10.pdf",height=6,width=10)
p5_ageComp <- ggboxplot(plot.info.age, 
                           x = "Gene_symbol", y = "Expression",
                           color = "black", fill = "Group",palette = "lancet",
                           xlab = "", ylab = "Gene expression",
                           main = "Gene expression comparison (age: <60 vs. >=60)")+ # age: <60 (203) vs. ??60 (296)
  
  stat_compare_means(aes(group = Group), label = "p.signif", method = "wilcox.test", hide.ns = F)+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30,hjust = 1))
  
  
  
  # stat_compare_means(label = "p.signif", method = "wilcox.test", # wilcox.test
  #                    ref.group = ".all.", hide.ns = F) +
  # theme_bw()+
  # theme(panel.background = element_blank(),
  #       panel.grid = element_blank(),
  #       axis.text.x = element_text(angle = 60,hjust = 1))
p5_ageComp

dev.off()
