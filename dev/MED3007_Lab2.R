###-------------------------------------------------------------------###
###     Data visualization and dimensional reduction with Rstudio     ###
###-------------------------------------------------------------------###

### SCHEDULE
###---------

### 1. Some useful function(s): apply (remember also table, summary)
### 2. Exercise
### 3. PCA: a small example
### BREAK (15 mins)
### 4. PCA: a larger genomics example
### 5. Exercise

## Some of these R exercises are adopted from "James et al. (2013).
## An Introduction to Statistical Learning with Applications in R"
## (Chapter 10: "Unsupervised Learning").

# Remember to specify YOUR OWN working directory here!
# setwd("~/Dropbox_UiO/Dropbox/MED3007_2023/day 3")
# ..sanity check: click on "Files" in the lower right corner window
## and check if you are in the correct place

###-----------------------###
### 1. USEFUL FUNCTION(S) ###
###-----------------------###

## the easiest way to do matrix & array operations is with the function
## "apply" (and all its variations)

## you can use the function "apply" to apply a function on columns or rows
## of a matrix without using an explicit loop. Syntax is as follows:
## apply(X, margin, function)
## "X" is the array (matrix),
## "margin" is the dimension that should be kept
## "function" is the function that is to be applied for each component

## remark: margin needs to be an available dimension!

# example
mat <- matrix(c(1:6), ncol=2, byrow=TRUE) # create a matrix
mat
apply(mat, 1, sum)
apply(mat, 2, sum)

## same commands exist for list objects: sapply/lapply
a1 <- 1:10
a2 <- 2:30
a3 <- 3:40
lst <- list(a1,a2,a3)
lst

lapply(lst, mean) # returns a list with the results
sapply(lst, mean) # returns an array/vector with the results


## what if we would like to apply some function stratified according to the
## levels of a factor? --> tapply

# tapply in a simple example
a1 <- 1:15
a1
myfac <- c(rep(1,5), rep(2,5), rep(3,5))
myfac
tapply(a1, myfac, mean) # returns the mean below each grouping variable (here the grouping variable is 1, 2 and 3)


# more complicated example: NCI60 data
library(ISLR)
nci.labs <- NCI60$labs # Sample labels (tissue type)
nci.data <- NCI60$data # Gene expression data set
# what if I would like to compute the mean of each gene within each tissue type?
tissue.means <- apply(nci.data, 2, function(x){tapply(x, nci.labs, mean)})
dim(tissue.means)
table(nci.labs)


###-------------###
### 2. EXERCISE ###
###-------------###

### Load the "pima.txt" data from day 2 of the course, and
### 1. use the function "apply" to compute the mean and variance of all variables
### 2. use the function "tapply" to compute the mean of the "BMI" for each
###    level of the variable "pregnant" (how many times the woman was pregnant)



###-------------------------###
### 3. PCA: A SMALL EXAMPLE ###
###-------------------------###

# Remember: the goal of PCA is to explain most of the variability in the data with a
# smaller number of variables than the original data set. It finds a low-dimensional
# representation of a data set that contains as much of the variation in the data as possible.

# The idea is that each of the n observations lives in p-dimensional space, however, not
# all of these dimensions are equally interesting. PCA tries to find the interesting ones.

## We start with a PCA on the USArrests data set, which is part of
## the base R package (see Lab 10.1 in James et al., 2013). We look at this data
## set first, because it is small (only has 4 variables) and is therefore
## convenient for understanding all the results of PCA.

# Let's have a look at the data set first.
?USArrests
states <- row.names(USArrests)
states

# some more statistics...
dim(USArrests)
names(USArrests)
head(USArrests)
summary(USArrests)
apply(USArrests, 2, mean)
apply(USArrests, MARGIN=2, FUN=var)

# Remember that the apply() function will apply a function (here mean() and var())
# to each row (MARGIN=1) or column (MARGIN=2) of a data set.


## We now perform principal components analysis using prcomp()
##------------------------------------------------------------

## By default, prcomp() centers the variables to have mean zero.

## By using the option scale=TRUE, we also scale the variables to have standard deviation 1.
## Why does it make sense here to scale the variables to have the same standard deviation?

pr.out <- prcomp(USArrests, scale=TRUE)

names(pr.out)
pr.out$center # variable's means, the ones used to center
pr.out$scale # variable's standard deviations, the ones used to scale
pr.out$sdev # standard deviations of the principal components (decreasing)
pr.out$rotation # loadings

pr.out$x # scores
dim(pr.out$x)

?prcomp

# proportion of explained variance (PEV) -> scree plot
dj <- pr.out$sdev
PEV <- dj^2 / sum(dj^2)
plot(PEV, type = 'b', xlab = 'number of components', ylab = 'PEV', main = 'USA arrests data')

## Seems like the first two components explain most of the variance

# interpret the components
barplot(pr.out$rotation[,1], ylim = c(-1, 1), main = 'First Component') # PC1
barplot(pr.out$rotation[,2], ylim = c(-1, 1), main = 'Second Component') # PC2
barplot(pr.out$rotation[,3], ylim = c(-1, 1), main = 'Third Component') # PC3

# Visualize PC1 and PC2
plot(pr.out$x[,1:2], xlab="PC 1", ylab=" PC 2", pch="")
text(pr.out$x[,1],pr.out$x[,2], states, cex=.7)

###-----------------------------------###
### 4. PCA: A LARGER GENOMICS EXAMPLE ###
###-----------------------------------###

## Now we move to a typical large-scale biological data set
## (this is partly based on Lab 10.3 in James et al., 2013).

## We have already seen the NCI60 cancer cell line microarray data set,
## consisting of 6830 gene expression measurements on 64 cancer cell lines.
## It is available in the R package ISLR, which is the compendium R package
## to the book by James et al. (2013).

# The same dataset as we already looked at briefly above
nci.labs <- NCI60$labs # Sample labels (tissue type)
nci.data <- NCI60$data # Gene expression data set

dim(nci.data)
table(nci.labs)

## NB: important to always check if the variables are on the columns before doing a PCA
## (they cannot be on the rows!) If they are *not* on the columns we can
## transpose the dataset with the function t(), i.e. you would run the code:
## nci.data <- t(nci.data)


## PCA on the NCI60 Data
##----------------------

## First, let us perform the PCA analysis after scaling the variables (genes)
## to have standard deviation one.

# PCA analysis after scaling the variables to standard deviation one:
pr.out <- prcomp(nci.data, scale=TRUE)

## Lastly, we calculate the proportion of variance explained (PVE), and visualise
## it via a scree plot. In addition, we also plot the cumulative proportion of variance
## explained cumsum (pve), which will reach 100% when all principal components are added up.

# Proportion of variance explained (PVE):
summary(pr.out)

# Calculate the proportion of variance explained (PVE) by hand,
# make a scree plot and plot the cumulative proportion of variance explained cumsum(pve):
pr.var <- pr.out$sdev^2
pve <- pr.var/sum(pr.var)
pve <- 100*pve

par(mfrow=c(1,2))
plot(pve,  type="o", ylab="PVE", xlab="Principal Component", col="blue")
plot(cumsum(pve), type="o", ylab="Cumulative PVE", xlab="Principal Component", col="brown3")

# How many principal components would you keep to achieve a good dimension reduction,
# while keeping most of the variability in the data set?

mysel80 <- which(cumsum(pve) > 80)[1] # explains 80% of the variability
mysel70 <- which(cumsum(pve) > 70)[1] # explains 70% of the variability

par(mfrow=c(1,2)) # plot contains two smaller plots next to each other
plot(pve,  type="o", ylab="PVE", xlab="Principal Component", col="blue")
abline(v = mysel80)
abline(v = mysel70, col=3)
plot(cumsum(pve), type="o", ylab="Cumulative PVE", xlab="Principal Component", col="brown3")
abline(v = mysel80)
abline(h = 80)
abline(v = mysel70, col=3)
abline(h = 70, col=3)

## If we decide to only keep the principal components that explains 70% of the variance,
## we end up with 24 components, which we can further analyse to better understand the
## relationships between the variables. For simplicity we only look at the first few components.
## We plot the first few principal component score vectors, to visualize the results.
## The observations (cell lines) corresponding to a given cancer type
## will be plotted in the same colour.

# we here define a "helper function", which assigns a different colour to each sample
# label (nci.labs) in the next plot
Cols=function(vec){
  cols=rainbow(length(unique(vec)))
  return(cols[as.numeric(as.factor(vec))])
}

# Plot the first vs second and first vs third principal component score vectors,
# with colors associated to labels (using the Cols() helper function)
par(mfrow=c(1,2))
plot(pr.out$x[,1:2], col=Cols(nci.labs), pch=19,xlab="PC 1",ylab=" PC 2")
plot(pr.out$x[,c(1,3)], col=Cols(nci.labs), pch=19,xlab="PC 1",ylab=" PC 3")
legend('topleft', col=rainbow(length(unique(nci.labs))), legend=unique(nci.labs), bty='n', lwd=2, cex=.6)


###-------------###
### 5. EXERCISE ###
###-------------###

## Consider again the gene expression data set "Ch10Ex11.csv"
## (which can be also found on the book website, www.StatLearning.com)
## that consists of 40 tissue samples with measurements on 1,000 genes.
## The first 20 samples are from healthy patients,
## while the second 20 are from a diseased group.

## 1. Load in the data using read.csv(). You will need to select header=F.
##    Alternatively: load in the data using "Import dataset" in the upper right window,
##    and click "no" on the "Heading" option.
## 2. Perform a PCA of these data and visualize the results. Note: remember to check
##    if the variables (genes) are on the columns in the dataset before running the PCA.
##    If they are not: use t() to transform the dataset.
