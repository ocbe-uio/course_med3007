###-------------###
### 2. EXERCISE ###
###-------------###

setwd("~/Dropbox_UiO/Dropbox/MED3007_2023/day 3")

### Load the "pima.txt" data from day 2 of the course, and
### 1. use the function "apply" to compute the mean and variance of all variables
### 2. use the function "tapply" to compute the mean of the "BMI" for each
###    level of the variable "pregnant" (how many times the woman was pregnant)


data <- read.table('./data/pima.txt', header = TRUE) # load in the data
head(data) # have a look at the dataset, does it look OK?

apply(data, 2, mean)
apply(data, 2, var)

tapply(data$bmi, data$pregnant, mean)


###-------------###
### 4. EXERCISE ###
###-------------###

## Consider again the gene expression data set "Ch10Ex11.csv"
## (which can be also found on the book website, www.StatLearning.com)
## that consists of 40 tissue samples with measurements on 1,000 genes.
## The first 20 samples are from healthy patients,
## while the second 20 are from a diseased group.

## 1. Load in the data using read.csv(). You will need to select header=F.
## 2. Perform a PCA of these data and visualize the results. Note: remember to check 
## if the variables (genes) are on the columns in the dataset, before running the PCA.


exp.data <-  read.csv("./data/Ch10Ex11.csv",header=FALSE)
# I want each row to represent a sample, and each column a gene
exp.data <- t(exp.data) 
dim(exp.data)
# should have n=40 samples/rows, and 1000 columns --> OK!
groups <- c(rep(1,20), rep(2,20)) # group variable

# PCA
pr.exp <- prcomp(exp.data, scale=TRUE)

# Plot proportion of variance explained
pr.var <- pr.exp$sdev^2
pve <- pr.var/sum(pr.var)
pve <- 100*pve
par(mfrow=c(1,2))
plot(pve, type="o", ylab="PVE", xlab="Principal Component", col="blue")
plot(cumsum(pve), type="o", ylab="Cumulative PVE", xlab="Principal Component", col="red")

## Looks like most of the principal components are needed to explain the data well. 
## Maybe we can decide to keep 25-30 components?

# Can also plot some of the first principal components

# Remember the use the helper-function to get colours
Cols=function(vec){
  cols=rainbow(length(unique(vec)))
  return(cols[as.numeric(as.factor(vec))])
}

par(mfrow=c(1,2)) # plot-window has two small plots next to each other
plot(pr.exp$x[,1:2], col=Cols(groups), pch=19, xlab="PC 1", ylab=" PC 2")
plot(pr.exp$x[,c(1,3)], col=Cols(groups), pch=19,xlab="PC 1",ylab=" PC 3")
legend('topleft', col=rainbow(length(unique(groups))), legend=paste('group ',unique(groups),sep=''), bty='n', lwd=2, cex=.6)

## The first component is evidently *very* important: this we see in both plots 
## (PVE and the PC1 vs PC2 and PC3). We also see that PC1 describes the groups.


