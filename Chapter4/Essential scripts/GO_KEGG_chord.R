setwd('/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/second paper/Pathway_analysis')
getwd()


#####################################################################
#######################3####### Data import #########################
#####################################################################
library(readxl)
EP_results<- read_excel("EP_results_rawp0.01Ndelt0.1_NEBC_mle_QN_RMA.xlsx")
View(EP_results)

library(tidyr)
EP_results<-EP_results %>% drop_na(`Gene name`)
EP_gene <- unique(EP_results[ , 2])
View(EP_gene)
openxlsx::write.xlsx(EP_gene,"parentGene_pathway.xlsx")
#####################################################################
############## Gene functional enrichment analysis ##################
#####################################################################
# Load database====
library(clusterProfiler)
library(org.Hs.eg.db)
database <- org.Hs.egSYMBOL2EG
database <- as.list(database)
# Conversion of gene IDs----
# Show entrez gene IDs
hub_genes = EP_gene
View(hub_genes)
nrow(hub_genes)
gene.list = hub_genes$`Gene name`
View(gene.list)
gene.list.enrich <- database[names(database) %in% gene.list]
View(gene.list.enrich)

# GO and KEGG enrichment====
edo.bp <- enrichGO(as.character(unlist(gene.list.enrich)), OrgDb = org.Hs.eg.db,
                   ont = 'BP', pAdjustMethod = 'BH', 
                   pvalueCutoff = 0.05, qvalueCutoff = 1,
                   keyType = 'ENTREZID')

edo.mf <- enrichGO(as.character(unlist(gene.list.enrich)), OrgDb = org.Hs.eg.db,
                   ont = 'MF', pAdjustMethod = 'BH', 
                   pvalueCutoff = 0.05, qvalueCutoff = 1,
                   keyType = 'ENTREZID')

edo.cc <- enrichGO(as.character(unlist(gene.list.enrich)), OrgDb = org.Hs.eg.db,
                   ont = 'CC', pAdjustMethod = 'BH', 
                   pvalueCutoff = 0.05, qvalueCutoff = 1,
                   keyType = 'ENTREZID')

edo.kegg <- enrichKEGG(as.character(unlist(gene.list.enrich)), organism = "hsa",
                       pAdjustMethod = 'BH', keyType="kegg",
                       pvalueCutoff = 0.05, qvalueCutoff = 1)


# gene enriched gene name (no need to modify for other cancers, this is just a developed function)----
id2gene <- function(id.list){
  tmp = bitr(strsplit(id.list,"/")[[1]],fromType = "ENTREZID",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
  return(paste(tmp$SYMBOL,collapse = "/"))
}

#GP&KEGG Results out====
## bp out
#edo.bp.out = as.data.frame(edo.bp)[1:100,]
edo.bp.out = as.data.frame(edo.bp)
edo.bp.out = subset(edo.bp.out, p.adjust < 0.05)
tmp = unlist(as.data.frame(sapply(edo.bp.out$geneID, id2gene))[1])
colnames(tmp) = NULL
edo.bp.out$geneName = tmp
openxlsx::write.xlsx(edo.bp.out,"AS.GO_BP.xlsx")

## cc out
#edo.cc.out = as.data.frame(edo.cc)[1:100,]
edo.cc.out = as.data.frame(edo.cc)
edo.cc.out = subset(edo.cc.out, p.adjust < 0.05)
tmp = unlist(as.data.frame(sapply(edo.cc.out$geneID, id2gene))[1])
colnames(tmp) = NULL
edo.cc.out$geneName = tmp
openxlsx::write.xlsx(edo.cc.out,"AS.GO_CC.xlsx")

## MF out
#edo.mf.out = as.data.frame(edo.mf)[1:100,]
edo.mf.out = as.data.frame(edo.mf)
edo.mf.out = subset(edo.mf.out, p.adjust < 0.05)
tmp = unlist(as.data.frame(sapply(edo.mf.out$geneID, id2gene))[1])
colnames(tmp) = NULL
edo.mf.out$geneName = tmp
openxlsx::write.xlsx(edo.mf.out,"AS.GO_MF.xlsx")

## KEGG out
#edo.kegg.out = as.data.frame(edo.kegg)[1:100,]
edo.kegg.out = as.data.frame(edo.kegg)
edo.kegg.out = subset(edo.kegg.out, p.adjust < 0.05)
tmp = unlist(as.data.frame(sapply(edo.kegg.out$geneID, id2gene))[1])
colnames(tmp) = NULL
edo.kegg.out$geneName = tmp
openxlsx::write.xlsx(edo.kegg.out,"AS.KEGG.xlsx")

###########################
##Bubble plot 
###########################

# Functions to draw plots (no need to modify for other cancers)----
DrawGOBubblePlot <- function(dat, category = "BP", top.number = 15, col="blue"){
  # Draw bubble plot using DAVID function enrichment results
  
  category = toupper(category)
  if (category == "BP"){
    main.title = "Biological Process"
  } else if (category == "CC"){
    main.title = "Cellular Components"
  } else if (category == "MF"){
    main.title = "Molecular Function"
  } else if (category == "KEGG"){
    main.title = "KEGG"
  } else {
    return("ERROR! Wrong input parameter [category].")
  }
  
  dat1 = dat[c(1:top.number),c(2,3,4,5,6,9)]
  dat1[,2] = str_remove(dat1[,2],"/.*")
  dat1[,3] = str_remove(dat1[,3],"/.*")
  dat1$Ratio = as.numeric(dat1$GeneRatio) / as.numeric(dat1$BgRatio)
  
  
  dat1$Description = capitalize(dat1$Description)
  dat1$Description = factor(dat1$Description,levels=dat1$Description[length(dat1$Description):1])
  dat1$pvalue = -log10(dat1$pvalue)# change here for displaying adjusted p or raw p
  dat1$adjpvalue = -log10(dat1$p.adjust)
  
  p = ggplot(dat1,aes(Ratio,Description)) +
    geom_point(aes(size=Count,colour=adjpvalue)) +
    scale_colour_gradient(low=col,high="red") + 
    labs(colour=expression(-log[10]("adj.P Value")),size="Gene counts",  
         x="Gene Ratio",y="",title=main.title) +
    theme_bw() +
    scale_x_continuous(limits = c(0,max(dat1$Ratio) * 1.2)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 50))
  
  return(p)
}
#====

# Read in data just saved above and generate the plots====

library(stringr)
library(utils)
library(Hmisc)
library(cowplot)
library(GOplot)
#BP
data = openxlsx::read.xlsx("AS.GO_BP.xlsx")
head(data)
p1 = DrawGOBubblePlot(data,"BP",10,"blue")
p1

#MF
data = openxlsx::read.xlsx("AS.GO_MF.xlsx")
head(data)
p2 = DrawGOBubblePlot(data,"MF",10,"blue")
p2

#CC
data = openxlsx::read.xlsx("AS.GO_CC.xlsx")
head(data)
p3 = DrawGOBubblePlot(data,"CC",10,"blue")
p3

#KEGG
data = openxlsx::read.xlsx("AS.KEGG.xlsx")
head(data)
p4 = DrawGOBubblePlot(data,"KEGG",10,"blue")
p4

## output plots together (Figure 2A-D)
plot.out = plot_grid(p1,p2,p3,p4,labels=c("A","B","C","D"),align = "h",nrow=2)
save_plot("AS_adjp_10.GO_KEGG.pdf",plot.out,base_width = 13, base_height = 9)
save_plot("AS_adjp_10.GO_KEGG.tiff",plot.out,base_width = 13, base_height = 9,
          compression = "lzw",dpi = 500)


#####################################################################
############################# Chord plot#############################
#####################################################################
# Chord chart
# https://www.jianshu.com/p/e4bb41865b7f?utm_campaign=hugo&utm_medium=reader_share&utm_content=note


#FOR 165 genes----
DEG_chord <- table_PC # DGE results came from DGE r file
dim(DEG_chord)
head(DEG_chord)

#import 165 hub gene GO BP results
ego_bp = openxlsx::read.xlsx("/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Gene identification/data_needed/165hubgenes.GO_BP.xlsx")
View(ego_bp)
class(ego_bp)
ego_bp_chord <- ego_bp[1:10,c(1,2,10,6)] # extract wanted columns
View(ego_bp_chord)
ego_bp_chord$geneName <- str_replace_all(ego_bp_chord$geneName,"/",",") 
names(ego_bp_chord)=c("ID","Term","Genes","adj_pval")
ego_bp_chord$Category = "BP"
head(ego_bp_chord)

genes_chord = data.frame(ID=DEG_chord$SYMBOL,
                         logFC=DEG_chord$logFC) # universal 22165 DEG； no need to subset for gene specifically in enrichment input(e.g., 'gene_PC_165genes')
head(genes_chord)
dim(genes_chord)

# #!!! Codes below are for check if genes_chord (22165)'s subset (for 165 genes as input for enrichGO) would get the same results/graph (YES, WILL!!!)
# genes_165_chord = subset(DEG_for_up_down_set_PC, DEG_for_up_down_set_PC$SYMBOL %in% topgenes_overlap_KME_tomIMC )
# genes_165_chord = data.frame(ID=genes_165_chord$SYMBOL,
#                          logFC=genes_165_chord$logFC) # 统一全部22165 DEG的就可以了； 不需要专门subset成做enrichment input的gene set (e.g., 'gene_PC_165genes')
# head(genes_165_chord)
# dim(genes_165_chord)

# Chord analysis setup
circ_165_go <- circle_dat(ego_bp_chord,genes_chord)
chord_165_go <- chord_dat(data=circ_165_go, genes=genes_chord,process = ego_bp_chord$Term) 
#GOChord(chord_165_go, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)
#GOChord(chord_165_go, space = 0.02, gene.order = 'logFC', gene.space = 0.2, gene.size = 4, process.label = 10, title = 'Gene and GO-term Network') # process.label = 10; 在方形窗口下legend能够在窗口内
Figure_165_go <- GOChord(chord_165_go, space = 0.02, gene.order = 'logFC', gene.space = 0.2, gene.size = 4, process.label = 10)+
  ggtitle('Gene and GO-term Network for 165 Genes in Green Module')+theme(plot.title = element_text(hjust = 0.5,colour = "black",face = "bold.italic",size = 20)) # process.label = 10; 在方形窗口下legend能够在窗口内

Figure_165_go # Figure 2E

pdf("Figure 2E.pdf", width = 12, height = 13)
Figure_165_go # Figure 2E
dev.off()





