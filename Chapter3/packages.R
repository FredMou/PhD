##################################################################
######################## package needed ##########################
##################################################################
#General Bioconductor packages
library(BiocManager)
library(Biobase)
library(oligoClasses)
library(affy)
#BiocManager::install("affxparser")
library(affxparser)
#Annotation and data import packages
library(ArrayExpress)
library(pd.clariom.d.human)
library(clariomdhumantranscriptcluster.db)
library(clariomdhumanprobeset.db)
#Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)
#Analysis and statistics packages
library(limma)
library(topGO)
library("reactome.db") # must install this BEFORE the "ReactomePA"
library(ReactomePA)
library(clusterProfiler)
#Plotting and color options packages
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)
#Formatting/documentation packages
library(rmarkdown)
library(BiocStyle)
library(dplyr)
library(tidyr)
#Helpers:
library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)
library(devtools)
library(reshape2)
library(FactoMineR)
library(factoextra) 
# For table format
library(kableExtra)
#===Libraries for WGCNA----
library(impute)
library(WGCNA);
options(stringsAsFactors = FALSE);
#installation codes:
#install.packages("xtable")
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
library(ggpubr)
library(readxl)
library(xlsx)
library(org.Hs.eg.db)
library(enrichplot)
library(stringr)
library(tibble)
library(ggupset)
library(DOSE)
#options(java.parameters = "-Xmx1000m") #https://code.google.com/archive/p/rexcel/issues/33
library(plyr)
library(gtools) # for fc->log2fc
library(ggpubr)
library(ggthemes)
library(venn)   
library(ggpolypath)
library(VennDiagram)
library(UpSetR)
library(GOplot)
library(ggthemes)
library(gmodels)
library(ggcorrplot)
library(circlize)
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)









#General Bioconductor packages
library(Biobase)
library(oligoClasses)
library(affy)
#BiocManager::install("affxparser")
library(affxparser)
#Annotation and data import packages
library(ArrayExpress)
library(pd.clariom.d.human)
library(clariomdhumantranscriptcluster.db)
library(clariomdhumanprobeset.db)

#Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)

#Analysis and statistics packages
library(limma)
library(topGO)
#library(ReactomePA)
library(clusterProfiler)

#Plotting and color options packages
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)

#Formatting/documentation packages
library(rmarkdown)
library(BiocStyle)
library(dplyr)
library(tidyr)

#Helpers:
library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)
library(devtools)

library(reshape2)
library(FactoMineR)
library(factoextra) 

# For table format
library(kableExtra)


#===Libraries for WGCNA----
library(WGCNA);
options(stringsAsFactors = FALSE);
#installation codes:
#install.packages("xtable")
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

library(ggpubr)
library(ggplot2)
library(clariomdhumantranscriptcluster.db)

###########################
##Bubble plot 
###########################
# Packages
library(Hmisc)
library(ggplot2)
library(stringr)
library(cowplot)











