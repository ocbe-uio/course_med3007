###--------------###
### 2. EXERCISE  ###
###--------------###

# Load the data from the file "Testfil_Rcourse.xlsx"
# and consider the variable "vitD_v1"
# 1. plot the histogram of this variable and save it, and then also plot the 
#    boxplot of vitD_v1 stratified according to gender.
#    Does this variable look normally distributed? test it.
# 2. If the data are gaussian, perform a t-test to verify
#    that vitD_v1 is different across gender groups (otherwise, which test?).

library(readxl)
Testfil_Rcourse <- read_excel("data/Testfil_Rcourse.xlsx")
View(Testfil_Rcourse)
Testfil_Rcourse <- as.data.frame(Testfil_Rcourse)

png('./Rplots/ex2_histogram_vitaminD.png', height = 500, width = 500)
hist(Testfil_Rcourse$vitD_v1, prob = TRUE, xlab = 'vitamin D', main = 'histogram of vitamin D')
dev.off() # for closing the plot file

png('./Rplots/ex2_boxplot_vitaminD.png', height = 500, width = 500)
boxplot(vitD_v1 ~ gender, Testfil_Rcourse)
dev.off() # for closing the plot file


vit1 <- Testfil_Rcourse[which(data$gender==0),'vitD_v1']
vit2 <- Testfil_Rcourse[which(data$gender==1),'vitD_v1']


qqnorm(vit1, main = 'vitamin D - females')
qqline(vit1)
qqnorm(vit2, main = 'vitamin D - males')
qqline(vit2)


shapiro.test(vit1)
shapiro.test(vit2) # data is normally distb.

t.test(vitD_v1 ~ gender, data) # there is a difference in vitD between groups (gender)





###-------------###
### 4. EXERCISE ###
###-------------###


## Consider the gene expression data set "Ch10Ex11.csv"
## that consists of 40 tissue samples with measurements on 1,000 genes.
## The first 20 samples are from healthy patients,
## while the second 20 are from a diseased group.

## 1. Load in the data using read.csv(). You will need to select header=F.

exp.data <-  read.csv("./data/Ch10Ex11.csv",header=FALSE)
# I want each row to represent a sample, and each column a gene
exp.data <- t(exp.data) 
dim(exp.data)
# should have n=40 samples/rows, and 1000 columns --> OK!

## 2. Have a look at the data and describe them with appropriate descriptive measures.
genes.means <- apply(exp.data, 2, mean)
genes.var <- apply(exp.data, 2, var)

plot(genes.means, xlab = 'genes', main = 'mean across samples', ylab = 'mean expr')
abline(h=0)

# we can try to compute the within group means and plot in different colors..
groups <- c(rep(1,20), rep(2,20))
genes.gr.means <- apply(exp.data, 2, function(x){tapply(x, groups, mean)})
plot(genes.gr.means[1,], xlab = 'genes', main = 'within-group mean across samples', 
     ylim = range(genes.gr.means), ylab = 'mean expr')
points(genes.gr.means[2,], col=2)
abline(h=0)

## 3. Your collaborator wants to know which genes differ the most across the two groups.
##    Suggest a way to answer this question, and apply it here.

## We have to apply t-test to all genes, and then correct for multiple testing

# let us first apply the t-test and get an idea of the significance in the data
# How many significant p-values (without any correction)
sum(pval.ttest < alpha)

# how many expected false positives?
Vexp <- alpha*dim(exp.data)[2]  # this is E(V) from the slides
Vexp

# Simple way of calculation adjusted p-values (using p.adjust()):
pval.fwer <- p.adjust(pval.ttest, method = "bonferroni")
pval.fdr <- p.adjust(pval.ttest, method = "BH")

# Number of significant p-values after Bonferroni correction
sum(pval.fwer < alpha) # conservative

# Number of significant p-values after BH correction
sum(pval.fdr  < alpha)

# Get the significant genes (here the genes are numbered from 1 to 1000)
sign.genes <- which(pval.fdr < alpha)

# Plot them
plot(genes.gr.means[1,], xlab = 'genes', main = 'significant genes in blue', 
     ylim = range(genes.gr.means), ylab = 'mean expr')
points(genes.gr.means[2,], col=2)
points(sign.genes, genes.gr.means[1,sign.genes], col=4, pch=4)
points(sign.genes,genes.gr.means[2,sign.genes], col=4, pch=4)
abline(h=0)
