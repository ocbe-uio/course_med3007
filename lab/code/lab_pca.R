# code for PCA


# NCI60  ----
## Now we move to a typical large-scale biological data set
## (this is partly based on Lab 10.3 in James et al., 2013).

## We have already seen the NCI60 cancer cell line microarray data set,
## consisting of 6830 gene expression measurements on 64 cancer cell lines.
## It is available in the R package ISLR, which is the compendium R package
## to the book by James et al. (2013).

library(ISLR)
nci.labs <- NCI60$labs # Sample labels (tissue type)
nci.data <- NCI60$data # Gene expression data set
# what if I would like to compute the mean of each gene within each tissue type?
tissue.means <- apply(nci.data, 2, function(x){tapply(x, nci.labs, mean)})
dim(tissue.means)
table(nci.labs)


## NB: important to always check if the variables are on the columns before doing a PCA
## (they cannot be on the rows!) If they are *not* on the columns we can
## transpose the dataset with the function t(), i.e. you would run the code:
## nci.data <- t(nci.data)


## PCA on the NCI60 Data

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





# Ch10Ex11 ----


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



