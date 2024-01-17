###----------------------------------------------###
###   MED3007: Summary of methods with RStudio   ###
###----------------------------------------------###

### OUTLINE
###---------

### 0. Statistical testing 
### 1. Multiple testing with correction (Lab 1)
### 2. Principal component analysis (Lab 2)
### 3. Clustering (Lab 3)
###   3.1. Hierarchical clustering
###   3.2. K-means clustering

# Set the working directory:
setwd("~/Dropbox_UiO/Dropbox/MED3007_2023/day 5")

###-------------------------###
### 0. STATISTICAL TESTING  ###
###-------------------------###

library(readxl) # reading from excel can be done using the "readxl" package, so we need to load it first
data_exc <- read_excel("data/Testfil_Rcourse.xlsx")
data_exc <- as.data.frame(data_exc) # remember to always force your excel-dataset to be a dataframe!

# A single t-test can be performed the following way (given that it is normally distributed)
out <- t.test(vitD_v1 ~ gender, data_exc) 
# ..where we test if vitD_v1 is different across gender groups
out$p.value # print out the p-value from the t-test
out$conf.int # can also get the confidence interval
out

# However, there are other tests you can do with similar syntax:
?wilcox.test() # nonparametric (non-normal data)
?chisq.test()
# .. and so on

# Summary statistics
summary(data_exc)


###----------------------###
### 1. MULTIPLE TESTING  ###
###----------------------###

## Steps:
## 1. Load data
## 2. Make sure you have the variables (that you want to test) column-wise (look at the data!!)
## 3. Compute the p-values for all the variables you want to test. Nice to visualize it with a histogram!
## 4. To do the correction:
##    a) Bonferroni: compute new (stricter) significance level: alpha/k, k=number of tests. 
##       You keep the original p-values < alpha/k (these are considered significant under Bonf. correction)
##    b) B-H: sort the p-values from small to large, compute the thresholds for each p-value. 
##       You keep the original p-values < threshold.

# Load in the data
data_csv <- read.csv("data/Ch10Ex11.csv", header=F) # data contains 40 samples of 1000 genes 

# Always good to have a look at the data
View(data_csv) 

# Remember to check if the variables (i.e genes) are on the columns! 
dim(data_csv) 

# Need to transform the data so we have the variables as columns
data_csv <- t(data_csv) 

# Grouping variable: the samples were grouped into healthy and diseased patients
groups <- c(rep(1,20), rep(2,20)) 

# Need to decide on a target significance level
alpha <- .05

# Apply t-test to all genes in each of the groups using "apply" and "t.test" -> histogram 
pval.ttest <- apply(data_csv, 2, function(x){t.test(x[which(groups==1)], x[which(groups==2)])$p.value})
hist(pval.ttest)

# Simple way of calculation adjusted p-values (using p.adjust()):
pval.fwer <- p.adjust(pval.ttest, method = "bonferroni")
pval.fdr <- p.adjust(pval.ttest, method = "BH")

# Number of significant p-values after Bonferroni correction
sum(pval.fwer < alpha) # conservative

# Number of significant p-values after BH correction
sum(pval.fdr  < alpha)

# Let us find the significant genes after correction
sign.genes.Bonf <- which(pval.fwer < alpha)
sign.genes.BH <- which(pval.fdr < alpha)

# Plot of the overall means with significant genes
genes.means <- apply(data_csv, 2, mean) # compute the means using "apply"
plot(genes.means, xlab = 'Genes', main = 'Mean across samples', ylab = 'Mean expression')
points(sign.genes.Bonf, genes.means[sign.genes.Bonf], col=3, pch=16) 
points(sign.genes.BH, genes.means[sign.genes.BH], col=4, pch=16)
legend('topright', col=3:4, bty='n', lwd=2, legend = c('Bonferroni','B-H')) # add colorcode 


###-------------------------###
### 2. PCA: GENOMIC EXAMPLE ###
###-------------------------###

## Steps:
## 1. Load data
## 2. Make sure you have the variables column-wise (look at the data!!)
## 3. Do PCA with "prcomp", and we typically set "scale=TRUE" inside "prcomp"
## 4. Visualize the results:
##    a) Compute and plot proportion of variance explained (PVE), try to decide on 
##       how many principal components you would choose. 
##    b) Plot the first principal components. What do you see?
##       Extra: add a color for the group labels (if they exist).

# Do PCA with "prcomp"
pr.out <- prcomp(data_csv, scale=TRUE)

# Proportion of variance explained
pr.var <- pr.out$sdev^2
pve <- pr.var/sum(pr.var)
pve <- 100*pve
par(mfrow=c(1,2))
plot(pve, type="o", ylab="PVE", xlab="Principal Component", col="blue")
plot(cumsum(pve), type="o", ylab="Cumulative PVE", xlab="Principal Component", col="red")

# How many principal components would you keep to achieve a good dimension reduction,
# while keeping most of the variability in the data set?
mysel80 <- which(cumsum(pve) > 80)[1] # explains 80% of the variability
mysel70 <- which(cumsum(pve) > 70)[1] # explains 70% of the variability

# Setting a threshold at 80% PVE is very common, which would here result in 29 
# components (much more managable than 1000!)
# Now what?
# If we decide to keep 29 components, we can for ex. detect highly expressed genes by 
# inspecting the different principal components (the very first principal components 
# should typically have high weights (both positive and negative weight) for the most 
# interesting variables/genes, and so on..). 

# Visualize results
# Helper function for colors (for the different groups)
Cols=function(vec){
  cols=rainbow(length(unique(vec)))
  return(cols[as.numeric(as.factor(vec))])
}

# Plot the different components (color for each group)
par(mfrow=c(1,2))
plot(pr.out$x[,1:2], col=Cols(groups), pch=19, xlab="PC 1", ylab=" PC 2")
plot(pr.out$x[,c(1,3)], col=Cols(groups), pch=19,xlab="PC 1",ylab=" PC 3")
legend('topleft', col=rainbow(length(unique(groups))), legend=paste('group ',unique(groups),sep=''), bty='n', lwd=2, cex=.6)


###---------------###
### 3. CLUSTERING ###
###---------------###

## Steps:
## 1. Load data
## 2. Make sure you have the variables column-wise (--> look at the data!!)
## 3.1 Hierarchical clustering:
##    a) Compute the distance matrix for your dataset using "dist" with a chosen distance measure, 
##        for ex. "euclidian" (but also good to try several distances)
##    b) Do the clustering with "hclust" with a linkage method, for ex. "complete" (but also here 
##        good to try several methods and chose the best at the end).
##    c) Plot the dendograms to decide which method to chose, and how many clusters to chose. 
##    d) Use "cutree" to cut the dendogram so that you finally get the clustering of the data.
##       Now you can f.ex compare the results from different methods, or, if you already have 
##       the group labels, you can now see if you were able to detect them!
## 3.2 K-means clustering:
##    a) Set a random seed 
##    b) Run clustering with "kmeans" where you need to specify: the number of clusters, and number 
##       of initial starts. 


## Hierarchical clustering
# We continue to use the same dataset, however with the samples shuffled 
# (since they are now ordered by group)
myshuffle <- sample(dim(data_csv)[1])
data_csv <- data_csv[myshuffle,]

# Calculate the distance matrix (default = Euclidean):
data.dist <- dist(data_csv)

# Alternatively, try all distances 
data.dist.e <- dist(data_csv, method="euclidean") # this is now exactly the same as data.dist
data.dist.c <- dist(data_csv, method="canberra")
data.dist.m <- dist(data_csv, method="manhattan")

# Perform clustering with different linkage methods:
hc.complete <- hclust(data.dist, method="complete")
hc.average <- hclust(data.dist, method="average")
hc.single <- hclust(data.dist, method="single")

## Compare the dendrograms and comment on the results in light of how
## the different linkage methods work.
par(mfrow=c(1,3))
plot(hc.single, labels=groups,  main="Single Linkage", xlab="", sub="")
plot(hc.complete, labels=groups, main="Complete Linkage", xlab="", sub="")
plot(hc.average, labels=groups, main="Average Linkage", xlab="", sub="")

# How to decide between linkage methods? We look for the most "compact" clusters.
par(mfrow=c(1,3))
plot(hc.single, main='Euclidian single', xlab='', labels=F,  sub='')
rect.hclust(hc.single, k=2)
rect.hclust(hc.single, k=3)
rect.hclust(hc.single, k=4)
plot(hc.complete, main='Euclidean complete',xlab='', labels=F, sub='')
rect.hclust(hc.complete, k=2)
rect.hclust(hc.complete, k=3)
rect.hclust(hc.complete, k=4)
plot(hc.average, main='Euclidean average',xlab='', labels=F, sub='')
rect.hclust(hc.average, k=2)
rect.hclust(hc.average, k=3)
rect.hclust(hc.average, k=4)

# We can see that we do not gain very much by having 3 clusters, compared to 2.. But we 
# can still compare the results.

# We use the function "cutree" to "cut" the dendrogram at a given number of groups
# Let's try with both k=2 and k=3 and see what is best
cluster.ec <- cutree(hc.complete, k=c(2,3)) # complete linkage, euclidian
cluster.ea <- cutree(hc.average, k=c(2,3))  # average linkage, euclidian
# Compare the two clustering results with each other with "table"
table(cluster.ec[,"2"], cluster.ea[,"2"]) # they completely agree!
table(cluster.ec[,"3"], cluster.ea[,"3"]) # here there are some minor disagreements --> conclusion: k=2!!

# How are the true labels distributed between clusters:
table(cluster.ec[,"2"], groups[myshuffle]) # completely correct, the groups are detected!

## K-means clustering
set.seed(4)
km.out2 <- kmeans(data_csv, 2, nstart=20) # 2 clusters
km.out3 <- kmeans(data_csv, 3, nstart=20) # 3 clusters

# We can print the cluster values
km.out2$cluster

# We can also directly compare the k-means result (along rows)
# with the hierarchical clustering result (along columns)
table(km.out2$cluster, cluster.ec[,"2"], deparse.level=2) # complete agreement!
table(km.out3$cluster, cluster.ec[,"3"], deparse.level=2)

# How are the true labels distributed between clusters:
table(km.out2$cluster, groups[myshuffle]) # completely correct, the groups are detected!

# Since there seems to be overall best agreement on the clusters when k=2 across all methods, we
# can we be fairly certain that this is the true clustering. However, in many other situations 
# you will not be so certain. So: you should try several methods, and chose the result that seems 
# to be most consistent across methods. 