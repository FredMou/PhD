# working directory set up
setwd('/Users/zhuofanmou/Documents/Phd/papers/paper2/OmicS_heatmap/DEAS_fr_EP')
getwd()




lname <- load("isoform_PSI_NEBC_mle_QN_RMA.rda")
lname

PSI_data_EP <- PSI
class(PSI_data_EP)
View(PSI_data_EP) 

library(readxl)
sigEvents <- read_excel("EP_results_rawp0.01Ndelt0.1_NEBC_mle_QN_RMA.xlsx")
View(sigEvents)

sigEvents <- dplyr::filter(sigEvents,  !is.na(`Gene name`)) # remove row with NA gene symbols
sigEvents <- sigEvents[order(sigEvents$`Splicing Pvalue`) , ]# Assending #sigEvents[order(-sigEvents$`Splicing Pvalue`) , ] # Descending


# top 20 based on p values without deirection of det(PSI) (i.e., including up&down regulated events)
top20_events <- sigEvents[c(1:20) , ]
View(top20_events)
write.table(top20_events,file="TOP20_sigAS_STATS.txt",sep="\t",col.names=F,quote=F)


top20_events_PSI <- subset(PSI_data_EP, rownames(PSI_data_EP) %in% top20_events$ProbeID) # == means match the order as well
View(top20_events_PSI)
DEAS_top20_PSI_forOmicShare <- rbind(ID=colnames(top20_events_PSI),top20_events_PSI)
View(DEAS_top20_PSI_forOmicShare)
write.table(DEAS_top20_PSI_forOmicShare,file="TOP20_sigAS_PSI.txt",sep="\t",col.names=F,quote=F)


# top 20 based on p values including direction of det(PSI) (i.e., top10 up and top10 down DEAS events)
## top 10 up
sigEvents_up <- sigEvents[sigEvents$`Delta PSI` > 0 , ]# up DEAS
View(sigEvents_up)
sigEvents_up <- sigEvents_up[order(sigEvents_up$`Splicing Pvalue`) , ]
top10_events_up <- sigEvents_up[c(1:10) , ]
View(top10_events_up)
write.table(top10_events_up,file="TOP10_up_sigAS_STATS.txt",sep="\t",col.names=F,quote=F)

top10_events_up_PSI <- subset(PSI_data_EP, rownames(PSI_data_EP) %in% top10_events_up$ProbeID) # == means match the order as well
View(top10_events_up_PSI)
DEAS_top10_up_PSI_forOmicShare <- rbind(ID=colnames(top10_events_up_PSI),top10_events_up_PSI)
View(DEAS_top10_up_PSI_forOmicShare)
write.table(DEAS_top10_up_PSI_forOmicShare,file="TOP10_up_sigAS_PSI.txt",sep="\t",col.names=F,quote=F)

## top 10 down
sigEvents_down <- sigEvents[sigEvents$`Delta PSI` < 0 , ]# down DEAS
View(sigEvents_down)
sigEvents_down <- sigEvents_down[order(sigEvents_down$`Splicing Pvalue`) , ]
top10_events_down <- sigEvents_down[c(1:10) , ]
View(top10_events_down)
write.table(top10_events_down,file="TOP10_down_sigAS_STATS.txt",sep="\t",col.names=F,quote=F)


top10_events_down_PSI <- subset(PSI_data_EP, rownames(PSI_data_EP) %in% top10_events_down$ProbeID) # == means match the order as well
View(top10_events_down_PSI)
DEAS_top10_down_PSI_forOmicShare <- rbind(ID=colnames(top10_events_down_PSI),top10_events_down_PSI)
View(DEAS_top10_down_PSI_forOmicShare)
write.table(DEAS_top10_down_PSI_forOmicShare,file="TOP10_down_sigAS_PSI.txt",sep="\t",col.names=F,quote=F)




#================================================= Boxplot for needed events





lname <- load("isoform_PSI_NEBC_mle_QN_RMA.rda")
lname

PSI_data_EP <- PSI
class(PSI_data_EP)
PSI_data_EP <- as.data.frame(PSI_data_EP)
View(PSI_data_EP) 

library(readxl)
sigEvents <- read_excel("EP_results_rawp0.01Ndelt0.1_NEBC_mle_QN_RMA.xlsx")
View(sigEvents)




## For ZWINT: TC1000010686.hg_1
ZWINT <- PSI_data_EP[rownames(PSI_data_EP) == 'TC1000010686.hg_1',]
View(ZWINT)
ZWINT <- as.data.frame(t(ZWINT))
ZWINT$id <- rownames(ZWINT)
ZWINT <- ZWINT[ c(2,1,4,3,6,5,8,7,10,9,11,12,13,14,15,16,17,18), ]
ZWINT <- ZWINT[ c(1,3,5,7,9,11,13,15,17,2,4,6,8,10,12,14,16,18), ]

ZWINT$Type <- c(rep('Benign', 9), rep('Malignant', 9))
ZWINT <- ZWINT[ , c(2,1,3)]


colnames(ZWINT)=c("id", "Expression", "Type")

rt <- ZWINT
#设置比较组
group=levels(factor(rt$Type))
rt$Type=factor(rt$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#每组(group)正态检验(Shapiro test)====
# https://www.datanovia.com/en/lessons/normality-test-in-r/
# #方法一
# Tu = rt[rt$Type=='<8',]
# N  = rt[rt$Type=='>=8',]
# shapiro_test(Tu$Expression)
# shapiro_test(N$Expression)

#方法二
rt %>% group_by(Type) %>% shapiro_test(Expression)
#====

#boxplot

## For doing test seperately
wilcox.test(rt$Expression[1:9], rt$Expression[10:18], paired = TRUE)
t.test(rt$Expression[1:9], rt$Expression[10:18], paired = TRUE)


# boxplot=ggboxplot(rt, x="Type", y="Expression", color="Type",
#                   xlab="",
#                   ylab="PSI",
#                   legend.title="Type",
#                
#                   add = "jitter",
#                   palette = c("#29377E", "#AA1F23"))+
#   
#   stat_compare_means(paired = TRUE, comparisons = my_comparisons,label = "p.signif", method = "t.test") # label = "p.signif" 加了就可以显示星星了； stat_compare_means() 默认是wilcox.test, 要改的话用“method = "t.test"”, 但结果差不多的
# 
# boxplot



boxplot=ggpaired(rt, x="Type", y="Expression", color="Type",
                  xlab="",
                  ylab="PSI",
                  legend.title="Type",
                  line.color = "grey", line.size = 0.4,
                  add = "jitter",
                  palette = c("#29377E", "#AA1F23"))+
  
  stat_compare_means(paired = TRUE, comparisons = my_comparisons,label = "p.signif", method = "t.test") # label = "p.signif" 加了就可以显示星星了； stat_compare_means() 默认是wilcox.test, 要改的话用“method = "t.test"”, 但结果差不多的

boxplot




#输出图片
pdf(file="ZWINT_inClariomD.pdf", width=5,height=4.5)
print(boxplot)
dev.off()




























