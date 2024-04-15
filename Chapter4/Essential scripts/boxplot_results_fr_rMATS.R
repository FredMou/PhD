library(ggpubr)
tumour <-c(0.333,0.259,0.077,0.189,0,0.429,0.101,0.091,0,0.081,0.129,0.259,0.278)
normal <-c(0.091,0.154,0,0,0.333,0,0.087,0.12,0.059,0.161,0,0.091,0.077)


######################### ZWINT Retained Intron, chr10:56358185-56358376", from rMATS
condition <- c("tumour","tumour","tumour","tumour","tumour","tumour","tumour","tumour","tumour","tumour","tumour","tumour","tumour",
               "normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal","normal")
PSI <- c(0.333,0.259,0.077,0.189,0,0.429,0.101,0.091,0,0.081,0.129,0.259,0.278,0.091,0.154,0,0,0.333,0,0.087,0.12,0.059,0.161,0,0.091,0.077)

d <- data.frame(condition = condition, PSI = PSI)
View(d)

ggpaired(d, x = "condition", y = "PSI", xlab = "Sample Status", ylab = "PSI value",
         color = "condition", line.color = "gray", line.size = 0.4,
         palette = "lancet",
         main = "ZWINT, Retained Intron, chr10:56358185-56358376, in PRJEB2449")+
  stat_compare_means(paired = TRUE, method = "wilcox.test") # wilcox.test




######################### ZWINT Retained Intron, chr10:56358185-56358376", from TCGA SpliceSeq
TCGAss_PSI <- load("PSI_TCGA.rda")
TCGAss_PSI
TCGAss_PSI <- PSIexp

View(TCGAss_PSI)
which(colnames(TCGAss_PSI) == "ZWINT|11811|RI")
ZWINT <- (TCGAss_PSI[ , c(15,16)])
View(ZWINT)

colnames(ZWINT)[2] <- "condition" 
class(ZWINT)
ZWINT <- as.data.frame(ZWINT)
ZWINT$condition[1:52] <- "normal"
ZWINT$condition[53:541] <- "tumour"

condition <- ZWINT$condition
length(condition)
PSI <- ZWINT$`ZWINT|11811|RI`
length(PSI)

d <- data.frame(condition = condition, PSI = PSI)
View(d)


p <- ggboxplot(d, x = "condition", y = "PSI", xlab = "Sample Status", ylab = "PSI value",
          color = "condition", palette = "lancet", add = "jitter")
p
p + stat_compare_means(method = "wilcox.test")




