###----------------------------------------------------------###
###     Data screening and multiple testing with RStudio     ###
###----------------------------------------------------------###

### OUTLINE
###---------
### 1. Recall of statistical testing in R (t-test) 
### 2. Exercise
### BREAK
### 3. Multiple testing correction
### 4. Exercise

setwd("C:/Users/valeriv/UiO Dropbox/Valeria Vitelli/Teaching/MED3007/MED3007_2024/day2")


###-----------------------------------------###
### 1. RECAP OF STATISTICAL TESTING: T-TEST ###
###-----------------------------------------###


## load and check the data first
load('./data/data_brainshake.RData')
data <- as.data.frame(mydatanew)
head(data)

# I use the data in the brainshake file
# Some descriptive statistics of each variable
summary(data)
names(data)

## descriptive plots: boxplot
# (for one specific variable across groups)
boxplot(BMI ~ kjonn, data)
boxplot(tot.kol ~ kjonn, data)
hist(data$BMI)

BMI_group1 <- data[which(data$kjonn==0),'BMI']
BMI_group2 <- data[which(data$kjonn==1),'BMI']

hist(BMI_group1) #men
hist(BMI_group2) #females

### 2 groups t-test
###---------------
# VERY STUPID TEST: is there a gender difference in BMI?

# t-test:
?t.test

out <- t.test(data[which(data$kjonn==0),'BMI'],data[which(data$kjonn==1),'BMI'], alternative='greater')
# much more clever syntax
out <- t.test(BMI ~ kjonn, data, alternative='greater')
out
# can also just run the t-test and not save it to an object
t.test(BMI ~ kjonn, data, alternative='greater')
# ..but often useful to save it!
# .. because then you can access the p-value and not have to copy-paste it from the console
names(out)
out$p.value
pvalue <- out$p.value
pvalue

## Here sample size is large (n = 162), so normality not needed
## BUT in general: how to check the normality assumption?

# QQ-plots:
par(mfrow=c(1,2)) # plot-window with two columns
qqnorm(BMI_group1, main = 'BMI - male')
qqline(BMI_group1)
qqnorm(BMI_group2, main = 'BMI - female')
qqline(BMI_group2)

# test data for normality
shapiro.test(BMI_group1)
shapiro.test(BMI_group2)

#### OTHER USEFUL TESTS: wilcox.test, prop.test, chisq.test

###-------------###
### 2. EXERCISE ###
###-------------###

# Load the data from the file "Testfil_Rcourse.xlsx"
# and consider the variable "vitD_v1"
# 1. plot the histogram of this variable and save it, and then also plot the 
#    boxplot of vitD_v1 stratified according to gender.
#    Does this variable look normally distributed?
# 2. Perform a t-test to verify that vitD_v1 is different across gender groups


###-------------------------------------------###
###    RECALL OF TESTING IN R - GENOMIC DATA  ###
###-------------------------------------------###

## the NCI60 data Example
##-----------------------

## Now we move to a typical large-scale biological data set
## which will be our running example throughout the course
## The NCI60 cancer cell line microarray data set consists of 6830 gene expression
## measurements on 64 cancer cell lines. It is available in the R package ISLR, 
## which is the compendium R package to the book by James et al. (2013).
## (it is not really a testing example, but we use it this way for illustration)

load('./data/NCI60.RData')

nci.labs <- NCI60$labs # Sample labels (tissue type)
nci.data <- NCI60$data # Gene expression data set
length(nci.labs)
dim(nci.data)
table(nci.labs)

# we have several cancer types but none of them has large dimension
# we thus create a grouping in 2 macro-groups of cancers, to test for 
# gene expression differences across these 2 groups only
ind1 <- which(nci.labs=='BREAST' | nci.labs=='OVARIAN' | nci.labs=='PROSTATE')
ind2 <- which(nci.labs=='LEUKEMIA' | nci.labs=='MELANOMA' | nci.labs=='NSCLC' | nci.labs=='RENAL') # blood, skin, lung, kidney
cancer.type <- c(rep(1, length(ind1)), rep(2, length(ind2))) # we have two "cancer type" groups
mydata <- nci.data[c(ind1, ind2) , ]
dim(mydata)
# those 2 groups should be relatively far in terms of genomic characteristics

# NB: make sure that the samples (patients) are on the rows, and variables (genes) are on the columns!


## First we look at the data and describe them with some descriptive measures
##---------------------------------------------------------------------------
# mean across patients for each gene
genes.means <- apply(mydata, 2, mean) #1=rows, 2=columns
genes.var <- apply(mydata, 2, var)

plot(genes.means, xlab = 'genes', main = 'mean across samples', ylab = 'mean expr')
abline(h=0, col=3)
# we can try to compute the within group means and plot in different colors..
genes.gr.means <- apply(mydata, 2, function(x){tapply(x, cancer.type, mean)})
plot(genes.gr.means[1,], xlab = 'genes', main = 'within-group mean across samples', 
     ylim = range(genes.gr.means), ylab = 'mean expr')
points(genes.gr.means[2,], col=2)
abline(h=0)


## Brute force application of the t-test across all genes
##-------------------------------------------------------
alpha <- .05
pval.ttest <- apply(mydata, 2, function(x){t.test(x[which(cancer.type==1)], x[which(cancer.type==2)])$p.value})

################
# Example: break down t-test for one gene:
x <- mydata[,1] #gene 1
x1 <- x[which(cancer.type==1)] #select values for group 1
x2 <- x[which(cancer.type==2)] #select values for group 2
res <- t.test(x1, x2) #do the t-test
res #show a summary of the test result
res$p.value #extract the p-value


hist(pval.ttest)

# bit more details for nicer plotting
hist(pval.ttest, breaks = 1/alpha, xlab = 'p-values', ylab = 'counts', 
     main = 'Histogram of t-test p-values across genes')
# the histogram of the p-values shows nearly uniform values
# but we see a slight peak around 0

# How many significant p-values (without any correction)
sum(pval.ttest < alpha)

# how many expected false positives?
Vexp <- alpha*dim(mydata)[2]  # this is the expected value of V from the slides
Vexp



###--------------------------------###
### 3. MULTIPLE TESTING CORRECTION ###
###--------------------------------###

# Simple way of calculation adjusted p-values (using p.adjust()):
pval.fwer <- p.adjust(pval.ttest, method = "bonferroni")
pval.fdr <- p.adjust(pval.ttest, method = "BH")

# Number of significant p-values after Bonferroni correction
sum(pval.fwer < alpha) # conservative

# Number of significant p-values after BH correction
sum(pval.fdr  < alpha)

# Conclusion: no significant difference in gene expression across the groups! 







###-------------###
### 4. EXERCISE ###
###-------------###

## Consider the gene expression data set "Ch10Ex11.csv"
## that consists of 40 tissue samples with measurements on 1,000 genes.
## The first 20 samples are from healthy patients,
## while the second 20 are from a diseased group.

## 1. Load in the data using read.csv(). You will need to select header=F.
## 2. Have a look at the data and describe them with appropriate descriptive measures.
## 3. Your collaborator wants to know which genes differ the most across the two groups.
##    Suggest a way to answer this question, and apply it here.




