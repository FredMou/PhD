
##################################################################
############################# data import ########################
##################################################################
# load finalised expression dataset from the pre-processing
lnames_limma_result <- load('/Users/zhuofanmou/All documents/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Gene identification/data_needed/table_PC.RData') # result output from limma
lnames_limma_result
DEGs_int <- read_excel("/Users/zhuofanmou/All documents/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Gene identification/data_needed/47_DGEonly_logtransformed_intensities.xlsx") # log2-transformed intensity for the 47 differentially expressed genes


##################################################################
###################### Histogram of raw p-value ##################
##################################################################
hist(table_PC$P.Value, col = brewer.pal(3, name = "Set2")[0],
     main = "Malignant vs Benign prostate tissues", xlab = "raw p-values")

# # For FDR p-value
# hist(table_PC$adj.P.Val, col = brewer.pal(3, name = "Set2")[0],
#      main = "Tumor vs Normal prostate tissues", xlab = "BH adjusted p-values")


##################################################################
############################# Volcano plot #######################
##################################################################
# Option1 (best)
DGE_results_PC = table_PC
DGE_results_PC$log_adj_p <- -log10(DGE_results_PC$adj.P.Val)
#adding new column 'group'
DGE_results_PC$Group = "not-significant"
DGE_results_PC$Group[which((DGE_results_PC$adj.P.Val<0.1) & (DGE_results_PC$logFC>log2(1.5)))] = "up-regulated"
DGE_results_PC$Group[which((DGE_results_PC$adj.P.Val<0.1) & (DGE_results_PC$logFC< -log2(1.5)))] = "down-regulated"
table(DGE_results_PC$Group)

DGE_results_PC$Label = ""
DGE_results_PC <- DGE_results_PC[order(DGE_results_PC$adj.P.Val),]
up.genes <- head(DGE_results_PC$SYMBOL[which(DGE_results_PC$Group == "up-regulated")], 10)
down.genes <- head(DGE_results_PC$SYMBOL[which(DGE_results_PC$Group == "down-regulated")], 10)
DEG_top10_PC <- c(as.character(up.genes), as.character(down.genes))
DGE_results_PC$Label[match(DEG_top10_PC, DGE_results_PC$SYMBOL)] <- DEG_top10_PC

# plot volcano
Figure_volc <- ggscatter(DGE_results_PC,
                         x = "logFC", 
                         y = "log_adj_p", 
                         color = "Group", 
                         palette = c("#2f5688", "#BBBBBB","#CC0000"), 
                         size = 1,
                         label = DGE_results_PC$Label,
                         font.label = 8,
                         repel = T,
                         xlab = "log2FoldChange",
                         ylab = "-log10(Adjust P-value)")+theme_base()+
  theme(plot.title = element_text(size=15, hjust = 0.5))+ # center)hjust = 0.5_ and change the font of title (size)
  ggtitle("Volcano plot, malignant vs. benign")+
  geom_hline(yintercept = -log10(0.1), #~'-log10(0.05)'
             linetype = "dashed")+
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), # log2(2) = 1
             linetype = "dashed")

Figure_volc




# Option2 (nice but not as good as option1)
DEG=table_PC
logFC_cutoff <- log2(1.5)  #LogFC Threshold
DEG$change = as.factor(ifelse(DEG$adj.P.Val < 0.1 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')    
                       #logFC>logFC_cutoff+adj.P<0.1->UP, if not then if abs(lofFC)>cutoff+adj.P<0.1->down, else->NOT;
)
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(DEG[DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(DEG[DEG$change =='DOWN',])
)
this_tile
#head(DEG)
g = ggplot(data=DEG, aes(x=logFC, y=-log10(adj.P.Val), color=change)) + 
  geom_point(alpha=0.4, size=1.75) + 
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + ylab("-log10 FDR corrected p-value") +
  ggtitle( this_tile  ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))  ## corresponding to the levels(res$change)
print(g)


##################################################################
################################ Heatmap #########################
##################################################################
DEGs_int <- as.data.frame(DEGs_int[ , c(2:20)])

rownames(DEGs_int) <- DEGs_int$SYMBOL 
DEGs_int <- DEGs_int[ , -19]
colnames(DEGs_int) <- c('10T','10B','12T','12B','13T','13B','14T','14B','15T','15B','17B','17T','7B','7T','8B','8T','9B','9T')

DEGs_int <- DEGs_int[ , c(1:10,12,11,14,13,16,15,18,17)]
View(DEGs_int)

Xtop <-as.matrix(DEGs_int)
View(Xtop)
colnames(Xtop) <- c("Patient10_Malignant",
                    "Patient10_Benign",
                    "Patient12_Malignant",
                    "Patient12_Benign",
                    "Patient13_Malignant",
                    "Patient13_Benign",
                    "Patient14_Malignant",
                    "Patient14_Benign",
                    "Patient15_Malignant",
                    "Patient15_Benign",
                    "Patient17_Malignant",
                    "Patient17_Benign",
                    "Patient7_Malignant",
                    "Patient7_Benign",
                    "Patient8_Malignant",
                    "Patient8_Benign",
                    "Patient9_Malignant",
                    "Patient9_Benign")

#heatmap(Xtop, Colv = NA, margins = c(7,2))
Heatmap_stadard <- pheatmap(Xtop, cluster_cols = F,cluster_rows = T, scale = "row", col = greenred(75), border_color = NA)
Heatmap_stadard


