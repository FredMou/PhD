
setwd("/Users/zhuofanmou/Documents/Phd/papers/paper2/Boxplot")
getwd()
#install.packages("ggpubr")

library(tidyverse)
library(ggpubr)
library(rstatix)

######################################################
############# t test was used for riskscore ##########
######################################################

#读取输入文件， 并对输入文件整理
# rt=read.table(inputFile, header=T, sep="\t", check.names=F)

# training ------------
lnames_train <- load('/Users/zhuofanmou/Library/CloudStorage/OneDrive-UniversityofExeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Model construction/R_data_needed/df_clinic.genemarkers.COX_train.rda')
lnames_train
rt <- df_clinic.genemarkers.COX_train
View(rt)


riskScore=colnames(rt)[8]
rt$id <- rownames(rt)
rt <- rt[ , c(20,8,19)]
colnames(rt)=c("id", "Expression", "Type")
rt[rt$Type == '<8',]$Type <- "<8"#"less than 8" # OR:0
rt[rt$Type == '≥8',]$Type <- ">=8"#"greater or equal to 8" # OR:1

#设置比较组
group=levels(factor(rt$Type))
rt$Type=factor(rt$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#每组(group)正态检验(Shapiro test)====
# https://www.datanovia.com/en/lessons/normality-test-in-r/
#方法一
Tu = rt[rt$Type=='<8',]
N  = rt[rt$Type=='>=8',]
shapiro_test(Tu$Expression)
shapiro_test(N$Expression)

#方法二
rt %>% group_by(Type) %>% shapiro_test(Expression)
#====

#boxplot
boxplot=ggboxplot(rt, x="Type", y="Expression", color="Type",
                  xlab="",
                  ylab="Five-gene risk score",
                  legend.title="Type",
                  
                  add = "jitter",
                  palette = c("#29377E", "#AA1F23"))+
  
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test") # label = "p.signif" 加了就可以显示星星了； stat_compare_means() 默认是wilcox.test, 要改的话用“method = "t.test"”, 但结果差不多的

boxplot

#输出图片
pdf(file="riskscore_Gleason_train.pdf", width=5,height=4.5)
print(boxplot)
dev.off()



#------------

rt <- df_clinic.genemarkers.COX_train
View(rt)


riskScore=colnames(rt)[8]
rt$id <- rownames(rt)
rt <- rt[ , c(20,8,17)]
colnames(rt)=c("id", "Expression", "Type")
rt <- rt %>% drop_na(Type)

#设置比较组
group=levels(factor(rt$Type))
rt$Type=factor(rt$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#每组(group)正态检验(Shapiro test)====
# https://www.datanovia.com/en/lessons/normality-test-in-r/
#方法一
Tu = rt[rt$Type=='T2',]
N  = rt[rt$Type=='T3-4',]
shapiro_test(Tu$Expression)
shapiro_test(N$Expression)

#方法二
rt %>% group_by(Type) %>% shapiro_test(Expression)
#====

#boxplot
boxplot=ggboxplot(rt, x="Type", y="Expression", color="Type",
                  xlab="",
                  ylab="Five-gene risk score",
                  legend.title="Type",
                  palette = c("#323232", "#989898"),#palette = c("#5588BA", "#E13446"),
                  add = "jitter")+
  
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test") # label = "p.signif" 加了就可以显示星星了； stat_compare_means() 默认是wilcox.test, 要改的话用“method = "t.test"”, 但结果差不多的

boxplot

#输出图片
pdf(file="riskscore_T_train.pdf", width=5,height=4.5)
print(boxplot)
dev.off()




#------------

rt <- df_clinic.genemarkers.COX_train
View(rt)


riskScore=colnames(rt)[8]
rt$id <- rownames(rt)
rt <- rt[ , c(20,8,18)]
colnames(rt)=c("id", "Expression", "Type")
rt <- rt %>% drop_na(Type)

#设置比较组
group=levels(factor(rt$Type))
rt$Type=factor(rt$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#每组(group)正态检验(Shapiro test)====
# https://www.datanovia.com/en/lessons/normality-test-in-r/
#方法一
Tu = rt[rt$Type=='N0',]
N  = rt[rt$Type=='N1',]
shapiro_test(Tu$Expression)
shapiro_test(N$Expression)

#方法二
rt %>% group_by(Type) %>% shapiro_test(Expression)
#====

#boxplot
boxplot=ggboxplot(rt, x="Type", y="Expression", color="Type",
                  xlab="",
                  ylab="Five-gene risk score",
                  legend.title="Type",
                  palette = c("#323232", "#989898"),#palette = c("#5588BA", "#E13446"),
                  add = "jitter")+
  
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test") # label = "p.signif" 加了就可以显示星星了； stat_compare_means() 默认是wilcox.test, 要改的话用“method = "t.test"”, 但结果差不多的

boxplot

#输出图片
pdf(file="riskscore_N_train.pdf", width=5,height=4.5)
print(boxplot)
dev.off()



#------------

rt <- df_clinic.genemarkers.COX_train
View(rt)


riskScore=colnames(rt)[8]
rt$id <- rownames(rt)
rt <- rt[ , c(20,8,16)]
colnames(rt)=c("id", "Expression", "Type")
rt <- rt %>% drop_na(Type)
rt[rt$Type == '<60',]$Type <- "<60"#"less than 8" # OR:0
rt[rt$Type == '≥60',]$Type <- ">=60"#"greater or equal to 8" # OR:1

#设置比较组
group=levels(factor(rt$Type))
rt$Type=factor(rt$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#每组(group)正态检验(Shapiro test)====
# https://www.datanovia.com/en/lessons/normality-test-in-r/
#方法一
Tu = rt[rt$Type=='<60',]
N  = rt[rt$Type=='>=60',]
shapiro_test(Tu$Expression)
shapiro_test(N$Expression)

#方法二
rt %>% group_by(Type) %>% shapiro_test(Expression)
#====

#boxplot
boxplot=ggboxplot(rt, x="Type", y="Expression", color="Type",
                  xlab="",
                  ylab="Five-gene risk score",
                  legend.title="Type",
                  palette = c("#323232", "#989898"),#palette = c("#5588BA", "#E13446"),
                  add = "jitter")+
  
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test") # label = "p.signif" 加了就可以显示星星了； stat_compare_means() 默认是wilcox.test, 要改的话用“method = "t.test"”, 但结果差不多的

boxplot

#输出图片
pdf(file="riskscore_age_train.pdf", width=5,height=4.5)
print(boxplot)
dev.off()


#------------

rt <- df_clinic.genemarkers.COX_train
View(rt)


riskScore=colnames(rt)[8]
rt$id <- rownames(rt)
rt <- rt[ , c(20,8,2)]
colnames(rt)=c("id", "Expression", "Type")
rt[rt$Type == 1,]$Type <- "Progressed"#"less than 8" # OR:0
rt[rt$Type == 0,]$Type <- "Not progressed"#"greater or equal to 8" # OR:1

#设置比较组
group=levels(factor(rt$Type))
rt$Type=factor(rt$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#每组(group)正态检验(Shapiro test)====
# https://www.datanovia.com/en/lessons/normality-test-in-r/
#方法一
Tu = rt[rt$Type=='Progressed',]
N  = rt[rt$Type=='Not progressed',]
shapiro_test(Tu$Expression)
shapiro_test(N$Expression)

#方法二
rt %>% group_by(Type) %>% shapiro_test(Expression)
#====






#boxplot
boxplot=ggboxplot(rt, x="Type", y="Expression", color="Type",
                  xlab="",
                  ylab="Five-gene risk score",
                  legend.title="Type",
                  palette = c("#323232", "#989898"),#palette = c("#5588BA", "#E13446"),
                  add = "jitter")+
  
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test") # label = "p.signif" 加了就可以显示星星了； stat_compare_means() 默认是wilcox.test, 要改的话用“method = "t.test"”, 但结果差不多的

boxplot

#输出图片
pdf(file="riskscore_progressionStatus_train.pdf", width=5,height=4.5)
print(boxplot)
dev.off()





# testing ------------
lnames_test <- load('/Users/zhuofanmou/OneDrive - University of Exeter/PhD work/my project/my papers/first paper/latest_docs/Paper1_draft_v2/more recent docs/R codes/Model construction/R_data_needed/df_clinic.genemarkers.COX_test.rda')
lnames_test
rt <- df_clinic.genemarkers.COX_test
View(rt)


riskScore=colnames(rt)[8]
rt$id <- rownames(rt)
rt <- rt[ , c(20,8,19)]
colnames(rt)=c("id", "Expression", "Type")
rt[rt$Type == '<8',]$Type <- "<8"#"less than 8" # OR:0
rt[rt$Type == '≥8',]$Type <- ">=8"#"greater or equal to 8" # OR:1

#设置比较组
group=levels(factor(rt$Type))
rt$Type=factor(rt$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#每组(group)正态检验(Shapiro test)====
# https://www.datanovia.com/en/lessons/normality-test-in-r/
#方法一
Tu = rt[rt$Type=='<8',]
N  = rt[rt$Type=='>=8',]
shapiro_test(Tu$Expression)
shapiro_test(N$Expression)

#方法二
rt %>% group_by(Type) %>% shapiro_test(Expression)
#====

#boxplot
boxplot=ggboxplot(rt, x="Type", y="Expression", color="Type",
                  xlab="",
                  ylab="Five-gene risk score",
                  legend.title="Type",
                  palette = c("#323232", "#989898"),#palette = c("#5588BA", "#E13446"),
                  add = "jitter")+
  
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test") # label = "p.signif" 加了就可以显示星星了； stat_compare_means() 默认是wilcox.test, 要改的话用“method = "t.test"”, 但结果差不多的

boxplot

#输出图片
pdf(file="riskscore_Gleason_test.pdf", width=5,height=4.5)
print(boxplot)
dev.off()



#------------

rt <- df_clinic.genemarkers.COX_test
View(rt)


riskScore=colnames(rt)[8]
rt$id <- rownames(rt)
rt <- rt[ , c(20,8,17)]
colnames(rt)=c("id", "Expression", "Type")
rt <- rt %>% drop_na(Type)

#设置比较组
group=levels(factor(rt$Type))
rt$Type=factor(rt$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#每组(group)正态检验(Shapiro test)====
# https://www.datanovia.com/en/lessons/normality-test-in-r/
#方法一
Tu = rt[rt$Type=='T2',]
N  = rt[rt$Type=='T3-4',]
shapiro_test(Tu$Expression)
shapiro_test(N$Expression)

#方法二
rt %>% group_by(Type) %>% shapiro_test(Expression)
#====

#boxplot
boxplot=ggboxplot(rt, x="Type", y="Expression", color="Type",
                  xlab="",
                  ylab="Five-gene risk score",
                  legend.title="Type",
                  palette = c("#323232", "#989898"),#palette = c("#5588BA", "#E13446"),
                  add = "jitter")+
  
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test") # label = "p.signif" 加了就可以显示星星了； stat_compare_means() 默认是wilcox.test, 要改的话用“method = "t.test"”, 但结果差不多的

boxplot

#输出图片
pdf(file="riskscore_T_test.pdf", width=5,height=4.5)
print(boxplot)
dev.off()




#------------

rt <- df_clinic.genemarkers.COX_test
View(rt)


riskScore=colnames(rt)[8]
rt$id <- rownames(rt)
rt <- rt[ , c(20,8,18)]
colnames(rt)=c("id", "Expression", "Type")
rt <- rt %>% drop_na(Type)

#设置比较组
group=levels(factor(rt$Type))
rt$Type=factor(rt$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#每组(group)正态检验(Shapiro test)====
# https://www.datanovia.com/en/lessons/normality-test-in-r/
#方法一
Tu = rt[rt$Type=='N0',]
N  = rt[rt$Type=='N1',]
shapiro_test(Tu$Expression)
shapiro_test(N$Expression)

#方法二
rt %>% group_by(Type) %>% shapiro_test(Expression)
#====

#boxplot
boxplot=ggboxplot(rt, x="Type", y="Expression", color="Type",
                  xlab="",
                  ylab="Five-gene risk score",
                  legend.title="Type",
                  palette = c("#323232", "#989898"),#palette = c("#5588BA", "#E13446"),
                  add = "jitter")+
  
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test") # label = "p.signif" 加了就可以显示星星了； stat_compare_means() 默认是wilcox.test, 要改的话用“method = "t.test"”, 但结果差不多的

boxplot

#输出图片
pdf(file="riskscore_N_test.pdf", width=5,height=4.5)
print(boxplot)
dev.off()



#------------

rt <- df_clinic.genemarkers.COX_test
View(rt)


riskScore=colnames(rt)[8]
rt$id <- rownames(rt)
rt <- rt[ , c(20,8,16)]
colnames(rt)=c("id", "Expression", "Type")
rt <- rt %>% drop_na(Type)
rt[rt$Type == '<60',]$Type <- "<60"#"less than 8" # OR:0
rt[rt$Type == '≥60',]$Type <- ">=60"#"greater or equal to 8" # OR:1

#设置比较组
group=levels(factor(rt$Type))
rt$Type=factor(rt$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#每组(group)正态检验(Shapiro test)====
# https://www.datanovia.com/en/lessons/normality-test-in-r/
#方法一
Tu = rt[rt$Type=='<60',]
N  = rt[rt$Type=='>=60',]
shapiro_test(Tu$Expression)
shapiro_test(N$Expression)

#方法二
rt %>% group_by(Type) %>% shapiro_test(Expression)
#====

#boxplot
boxplot=ggboxplot(rt, x="Type", y="Expression", color="Type",
                  xlab="",
                  ylab="Five-gene risk score",
                  legend.title="Type",
                  palette = c("#323232", "#989898"),#palette = c("#5588BA", "#E13446"),
                  add = "jitter")+
  
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test") # label = "p.signif" 加了就可以显示星星了； stat_compare_means() 默认是wilcox.test, 要改的话用“method = "t.test"”, 但结果差不多的

boxplot

#输出图片
pdf(file="riskscore_age_test.pdf", width=5,height=4.5)
print(boxplot)
dev.off()


#------------

rt <- df_clinic.genemarkers.COX_test
View(rt)


riskScore=colnames(rt)[8]
rt$id <- rownames(rt)
rt <- rt[ , c(20,8,2)]
colnames(rt)=c("id", "Expression", "Type")
rt[rt$Type == 1,]$Type <- "Progressed"#"less than 8" # OR:0
rt[rt$Type == 0,]$Type <- "Not progressed"#"greater or equal to 8" # OR:1

#设置比较组
group=levels(factor(rt$Type))
rt$Type=factor(rt$Type, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#每组(group)正态检验(Shapiro test)====
# https://www.datanovia.com/en/lessons/normality-test-in-r/
#方法一
Tu = rt[rt$Type=='Progressed',]
N  = rt[rt$Type=='Not progressed',]
shapiro_test(Tu$Expression)
shapiro_test(N$Expression)

#方法二
rt %>% group_by(Type) %>% shapiro_test(Expression)
#====






#boxplot
boxplot=ggboxplot(rt, x="Type", y="Expression", color="Type",
                  xlab="",
                  ylab="Five-gene risk score",
                  legend.title="Type",
                  palette = c("#323232", "#989898"),#palette = c("#5588BA", "#E13446"),
                  add = "jitter")+
  
  stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test") # label = "p.signif" 加了就可以显示星星了； stat_compare_means() 默认是wilcox.test, 要改的话用“method = "t.test"”, 但结果差不多的

boxplot

#输出图片
pdf(file="riskscore_progressionStatus_test.pdf", width=5,height=4.5)
print(boxplot)
dev.off()
