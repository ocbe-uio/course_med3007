###----------------------------------------------------------###
###     Data screening and multiple testing with RStudio     ###
###----------------------------------------------------------###


# setwd("~/Dropbox_UiO/Dropbox/MED3007_2023/day 2")


### descriptive statistics plots
###-----------------------------

## load and check the data first
load('./data/data_brainshake.RData')
data <- as.data.frame(mydatanew)
head(data)

# I use the data in the brainshake file
# descriptive statistics of each variable
summary(data)
names(data)

## descriptive plots: boxplot
# (for one specific variable across groups)
boxplot(BMI ~ kjonn, data)
boxplot(tot.kol ~ kjonn, data)
boxplot(triglyserider ~ kjonn, data)

par(mfrow=c(1,3))
boxplot(triglyserider ~ kjonn, data, subset = which(time=='t1'),main='First time point')
boxplot(triglyserider ~ kjonn, data, subset = which(time=='t2'),main='Second time point')
boxplot(triglyserider ~ kjonn, data, subset = which(time=='t3'),main='Third time point')

### 2 groups t-test
###---------------
# VERY STUPID TEST: is there a gender difference in BMI?

# t-test:
?t.test

out <- t.test(data[which(data$kjonn==0),'BMI'],data[which(data$kjonn==1),'BMI'], alternative='greater')
# much more clever syntax
out <- t.test(BMI ~ kjonn, data, alternative='greater')
names(out)
pvalue <- out$p.value
pvalue
# make the plot nicer and put the test result
par(mfrow=c(1,1))
boxplot(BMI ~ kjonn, data, axes = FALSE, ylab='BMI',
        main = 'BMI comparison between genders', ylim=range(data$BMI))
axis(1, labels=c('M','F'), at=1:2)
axis(2)
segments(1,max(data$BMI),2,max(data$BMI),lwd=2)
segments(1,max(data$BMI)-1,1,max(data$BMI),lwd=2)
segments(2,max(data$BMI)-1,2,max(data$BMI),lwd=2)
text(1.5,max(data$BMI)+.3,paste('p-value = ',round(pvalue,4),sep=''))


### PAIRED 1 group t-test
###----------------------
# This applies when the grouping is given by time (same sample across time points)
par(mfrow=c(1,1))
boxplot(X18_3 ~ time, data)

# I select the first and last time point
out <- t.test(data[which(data$time=='t1'),'X18_3'],
              data[which(data$time=='t3'),'X18_3'], alternative='greater',
              paired = TRUE)
out

#### SEE ALSO (USEFUL TESTS): prop.test, chisq.test


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

BiocManager::install("ISLR", update=FALSE)
library(ISLR)

# if this does not work:
# load('./data/NCI60.RData')

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


## First we look at the data and describe them with some descriptive measures
##---------------------------------------------------------------------------
genes.means <- apply(mydata, 2, mean) #1=rows, 2=columns
genes.var <- apply(mydata, 2, var)

plot(genes.means, xlab = 'genes', main = 'mean across samples', ylab = 'mean expr')
abline(h=0, col=3)
# interesting!
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

#Manual way of doing the multiple testing adjustment ourselves:

# how many rejections?
sum(pval.ttest < alpha) # number of rejections (R in the slides of Lecture 1)
# how many expected false positives?
Vexp <- alpha*dim(mydata)[2]  # this is E(V) from the slides
Vexp

## Estimate of the false discovery proportion
## (see also 6.9.1 of	Holmes & Huber (2019))
pi0 <- 2 * mean(pval.ttest > 0.5)
hist(pval.ttest, breaks = 1/alpha, xlab = 'p-values', ylab = 'counts',
     main = 'Histogram of t-test p-values across genes')
abline(h = pi0 * Vexp, lwd=2, col='blue')
abline(v = alpha, lwd=2, col='red')
# all the rejections that lie to the left of the red line and above the blue line
# are all significant hypotheses. How many are they?
sum(pval.ttest < alpha) - pi0 * Vexp
# but we cannot in this way understand WHICH ARE the significant ones..

###--------------------------------###
### 3. MULTIPLE TESTING CORRECTION ###
###--------------------------------###

## decide on a target significance level
alpha <- .05
k <- dim(mydata)[2]

## Bonferroni correction
adj.alpha <- alpha/k
# number of rejections when controlling for FWER < alpha
sum(pval.ttest < adj.alpha)
# Bonferroni is very conservative!
min(pval.ttest)
adj.alpha

## Benjamini-hockberg procedure
## (as described in slide 21, lecture 1;
##  also well described in 6.9.2 of Holmes & Huber (2019))

# first sort the p-values
pval.sorted <- sort.int(pval.ttest, index.return = TRUE, decreasing = FALSE)
# then compute the thresholds
b.h.thresh <- alpha*(1:k)/k
sum(pval.sorted$x <= b.h.thresh)

## Understanding a bit the B-H procedure
nplot <- 30 # I plot only the first 30 p-values
plot(1:nplot, b.h.thresh[1:nplot], type = 'l', ylim=range(pval.sorted$x[1:nplot]),
     xlab = 'ordered hypotheses', ylab = 'p-values',
     main = 'Benjamini & Hochberg adjustment')
points(1:nplot, pval.sorted$x[1:nplot], col=2, pch=16)
legend('topleft', col=1:2, bty = 'n', lwd=2,
       legend = c('B-H thresholds','ordered p-avalues'))

## which are the significant p-values? (in case there are any....)
# Bonferroni
sign.genes.Bonf <- which(pval.ttest < adj.alpha)
# B-H correction
sign.genes.BH <- pval.sorted$ix[which(pval.sorted$x <= b.h.thresh)]

# plot of the overall means with significant genes
plot(genes.means, xlab = 'genes', main = 'mean across samples', ylab = 'mean expr')
points(sign.genes.Bonf, genes.means[sign.genes.Bonf], col=3, pch=16)
points(sign.genes.BH, genes.means[sign.genes.BH], col=4, pch=16)
legend('topleft', col=3:4, bty='n', lwd=2, legend = c('Bonferroni','B-H'))
abline(h=0)

# nothing is significant...

###----------------------###
### 4. PERMUTATION TESTS ###
###----------------------###

## we have 47 samples, round 20 for each group... the sample size is
## quite low, and the t-test is maybe not justified.
## Let us try using a permutation test instead!

# decide on the number of permutations & significance level
alpha <- .05
B <- 100

# first compute the original t statistic
stat.ttest <- apply(mydata, 2, function(x){t.test(x[which(cancer.type==1)], x[which(cancer.type==2)])$statistic})
# then loop the permutation procedure
# (as described in slide 9, lecture 1)
stat.mat <- matrix(NA, B, dim(mydata)[2])
for(b in 1:B){
  # randomize the group assignment
  new.cancer.type <- rep(2, length(cancer.type))
  new.cancer.type[sample.int(length(cancer.type), length(ind1))] <- 1
  # compute permuted groups test statistic
  stat.ttest.perm <- apply(mydata, 2, function(x){t.test(x[which(new.cancer.type==1)], x[which(new.cancer.type==2)])$statistic})
  # store it in the matrix
  stat.mat[b, ] <- stat.ttest.perm
}

## for each gene, check whether the test statistic obtained with the original genes
## was EXTREME when compared to these test statistics. How?
# first compute the matrix absolute value
stat.abs.mat <- abs(stat.mat)
# then compare all 100 permutations with the original one,
# and count how many are larger (formula in the bottom of slide 9)
perm.p.val <- colSums(scale(stat.abs.mat, center = abs(stat.ttest), scale = FALSE) > 0 )/B

# let us compute R (number of rejections)
sum(perm.p.val < alpha)

# and the number of expected false positives E(V)
Vexp <- alpha*length(perm.p.val)
Vexp

# same visualization as before, but with the new p-values
pi0 <- 2 * mean(perm.p.val > 0.5)
hist(perm.p.val, breaks = 1/alpha, xlab = 'p-values', ylab = 'counts',
     main = 'Histogram of PERMUTED t-test p-values across genes')
abline(h = pi0 * Vexp, lwd=2, col='blue')
abline(v = alpha, lwd=2, col='red')

# How many should be the significant rejections now (on average)?
sum(perm.p.val < alpha) - pi0 * Vexp

# do we manage to reject something in this case? Bonferroni correction
adj.alpha <- alpha/k
sum(perm.p.val < adj.alpha)
# now we reject less, but stronger rejections!!
min(perm.p.val)
adj.alpha

# Benjamini-hockberg procedure
pval.sorted <- sort.int(perm.p.val, index.return = TRUE, decreasing = FALSE)
# then compute the thresholds
b.h.thresh <- alpha*(1:k)/k
sum(pval.sorted$x <= b.h.thresh)



## which are the significant p-values? (in case there are any....)
# Bonferroni
sign.genes.Bonf <- which(perm.p.val < adj.alpha)
# B-H correction
sign.genes.BH <- pval.sorted$ix[which(pval.sorted$x <= b.h.thresh)]

# plot of the overall means with significant genes
plot(genes.means, xlab = 'genes', main = 'mean across samples', ylab = 'mean expr')
points(sign.genes.Bonf, genes.means[sign.genes.Bonf], col=3, pch=16)
points(sign.genes.BH, genes.means[sign.genes.BH], col=4, pch=16)
legend('topleft', col=3:4, bty='n', lwd=2, legend = c('Bonferroni','B-H'))
abline(h=0)







