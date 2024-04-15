setwd('/Users/zhuofanmou/Documents/PhD Exeter/ClariomD/alternative_splicing/R_workspace_v2/survival_analysis')
getwd()

#==== PCaDB
PCaDB_clinical <- read.delim("~/Documents/PhD Exeter/ClariomD/alternative_splicing/R_workspace_v2/survival_analysis/PCaDB_clinical.txt")
View(PCaDB_clinical)
class(PCaDB_clinical)
PCaDB_clinical_tumour <- PCaDB_clinical[PCaDB_clinical$sample_type=='Primary' , ]
rownames(PCaDB_clinical_tumour) <- PCaDB_clinical_tumour$patient_id
View(PCaDB_clinical_tumour)

#==== TCGAbiolinks (limour)
TCGA.PRAD_clinical <- read.csv("~/Documents/PhD Exeter/ClariomD/alternative_splicing/R_workspace_v2/survival_analysis/TCGA-PRAD_clinical.csv")
rownames(TCGA.PRAD_clinical) <- TCGA.PRAD_clinical$bcr_patient_barcode
class(TCGA.PRAD_clinical)
View(TCGA.PRAD_clinical)

#==== PFS from Paper1
PFS <- load('TCGA_PAN_PFS_GSig_Input_edgeR.rda')
PFS <- rt_final_PFS_edgeR
View(PFS)
class(PFS)
table(PFS$vital_status)
PFS$ID <- rownames(PFS)
PFS$ID<-gsub("-01","",PFS$ID)
rownames(PFS) <- PFS$ID


#==== PFS and other clinicao-patho variables
clinical_final_PFS_forAS <- read.delim("clinical_final_PFS_forAS.txt")
View(clinical_final_PFS_forAS)
clinical_final_PFS_forAS <- clinical_final_PFS_forAS[ , -c(6,7)]


#==== TCGAbiolinks (limour) and other clinico-patho variables
TCGA.PRAD_clinical_modified <- read.delim("TCGA-PRAD_clinical_modified.txt")
View(TCGA.PRAD_clinical_modified)
rownames(TCGA.PRAD_clinical_modified) <- TCGA.PRAD_clinical_modified$bcr_patient_barcode
TCGA.PRAD_clinical_modified <- TCGA.PRAD_clinical_modified[ , c(171:187)]



#Combine1 (most important)
same=intersect(row.names(clinical_final_PFS_forAS),row.names(TCGA.PRAD_clinical_modified)) 
length(same)


d1=clinical_final_PFS_forAS[same,] 
dim(d1)
head(rownames(d1))
d2 <- TCGA.PRAD_clinical_modified[same,] 
dim(d2)
head(rownames(d2))

##Check if row names matched
table(row.names(d1) == row.names(d2))

d3 <- cbind(d1,d2) 
combined <- d3
# combined$patient_id <- rownames(combined)
View(combined)
write.table(combined,"PFSpaper1vsTCGAbiolinkslimour_allClinicopathoVars.txt",quote=F,row.names = T,col.names = T,sep="\t")

#Combine2
same=intersect(row.names(PCaDB_clinical_tumour),row.names(TCGA.PRAD_clinical)) 
length(same)


d1=PCaDB_clinical_tumour[same,] 
dim(d1)
head(rownames(d1))
d2 <- TCGA.PRAD_clinical[same,] 
dim(d2)
head(rownames(d2))

##Check if row names matched
table(row.names(d1) == row.names(d2))

d3 <- cbind(d1,d2) 
combined <- d3
# combined$patient_id <- rownames(combined)
View(combined)
write.table(combined,"PCaDBvsTCGAbiolinkslimour.txt",quote=F,row.names = F,col.names = T,sep="\t")



#Combine3
same=intersect(row.names(PFS),row.names(TCGA.PRAD_clinical)) 
length(same)


d1=PFS[same,] 
dim(d1)
head(rownames(d1))
d2 <- TCGA.PRAD_clinical[same,] 
dim(d2)
head(rownames(d2))

##Check if row names matched
table(row.names(d1) == row.names(d2))

d3 <- cbind(d1,d2) 
combined <- d3
# combined$patient_id <- rownames(combined)
View(combined)
write.table(combined,"PFSpaper1vsTCGAbiolinkslimour.txt",quote=F,row.names = F,col.names = T,sep="\t")


# Cmobine4
same=intersect(row.names(PFS),row.names(PCaDB_clinical_tumour)) 
length(same)


d1=PFS[same,] 
dim(d1)
head(rownames(d1))
d2 <- PCaDB_clinical_tumour[same,] 
dim(d2)
head(rownames(d2))

##Check if row names matched
table(row.names(d1) == row.names(d2))

d3 <- cbind(d1,d2) 
combined <- d3
# combined$patient_id <- rownames(combined)
View(combined)
write.table(combined,"PCaDBvsPFSpaper1.txt",quote=F,row.names = F,col.names = T,sep="\t")













