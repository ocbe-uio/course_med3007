###------------------------------###
###   MED3007: Exam simulation   ###
###------------------------------###

## Load exam data set using "load()"
load("data/MED3007_exam_data_V21.Rdata")

dim(clin)
dim(expr)
head(clin)

## Remember: we want to have the variables on the columns in our datasets, i.e. we want
## there to be 72 rows (the number of patients) for both datasets "expr" and "clin". 
## For "expr" we should have 7129 columns, since they are the genes, a.k.a our variables. 

# Use function t() to transpose the dataset so that it is on the correct form
data <- t(expr) # rename "expr" to "data"

# Make group variable (here that is whether a patient has ALL or AML)
groups <- clin$ALL.AML

## There are 2 groups: one possible analysis -> testing
##-----------------------------------------------------

# Apply t-test to all genes in each of the groups using "apply" and "t.test" -> histogram 
alpha <- 0.05
pval.ttest <- apply(data, 2, function(x){t.test(x[which(groups=="ALL")], x[which(groups=="AML")])$p.value})
hist(pval.ttest)

# Simple way of calculation adjusted p-values (using p.adjust()):
pval.fwer <- p.adjust(pval.ttest, method = "bonferroni")
pval.fdr <- p.adjust(pval.ttest, method = "BH")

# Number of significant p-values after Bonferroni correction
sum(pval.fwer < alpha) # conservative

# Number of significant p-values after BH correction
sum(pval.fdr  < alpha) # less conservative

# Let us find the significant genes after correction
sign.genes.Bonf <- which(pval.fwer < alpha)
sign.genes.BH <- which(pval.fdr < alpha)

# Plot of the overall means with significant genes
genes.means <- apply(data, 2, mean) # compute the means using "apply"
plot(genes.means, xlab = 'Genes', main = 'Mean across samples', ylab = 'Mean expression')
points(sign.genes.Bonf, genes.means[sign.genes.Bonf], col=3, pch=16) 
points(sign.genes.BH, genes.means[sign.genes.BH], col=4, pch=16)
legend('bottomright', col=3:4, bty='n', lwd=2, legend = c('Bonferroni','B-H')) # add colorcode 




## There are 2 groups: another possible analysis -> clustering
##--------------------------------------------------------

# Calculate the distance matrix, here with the Euclidean distance
data.dist <- dist(data, method = "euclidean")

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

# How to decide between linkage methods? We look for the most "compact" clusters, where 
# the division between the new branches are a bit "balanced".
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

## The different linkages seems to give very different results!
## Complete linkage seems to give 2 clusters/groups, which we know to be the truth. 
## It is also reassuring to see that when we try to divide the patients into more than 2 groups, 
## we do not gain anything. Although the other two linkage methods seems to say there 
## are no groups at all, only one big group. 
## Let's see if we detect the true groups for the clustering result from complete linkage. 

# We use the function "cutree" to "cut" the dendrogram into both k=2 and k=3 groups
cluster.ec <- cutree(hc.complete, k=c(2,3)) # complete linkage, euclidian

# How are the true labels distributed between clusters:
table(cluster.ec[,"2"], groups) # look only at k=2

# Hierarchical clustering was not able to detect the groups very well.. Let's try K-Means.

## K-means clustering
set.seed(4)
km.out2 <- kmeans(data, 2, nstart=20) # 2 clusters
km.out3 <- kmeans(data, 3, nstart=20) # 3 clusters

# We can directly compare the k-means result (along rows)
# with the hierarchical clustering result (along columns)
table(km.out2$cluster, cluster.ec[,"2"], deparse.level=2) 
table(km.out3$cluster, cluster.ec[,"3"], deparse.level=2)
# Not too bad, but still some disagreement

# How are the true labels distributed between clusters:
table(km.out2$cluster, groups) 

# Not very good here either, maybe even worse than hierarchcial clustering..?


## Also dimensional reduction via PCA can be an option
##----------------------------------------------------

pr.out <- prcomp(data, scale=TRUE)

# Plot also the proportion of variance explained
# Proportion of variance explained
pr.var <- pr.out$sdev^2
pve <- pr.var/sum(pr.var)
pve <- 100*pve
par(mfrow=c(1,2))
plot(pve, type="o", ylab="PVE", xlab="Principal Component", col="blue")
plot(cumsum(pve), type="o", ylab="Cumulative PVE", xlab="Principal Component", col="red")

# Can also try to plot only the first 10 components in the first plot to zoom in a bit
par(mfrow=c(1,1))
plot(1:10,pve[1:10], type="o", ylab="PVE", xlab="Principal Component", col="blue")
# -> the first three components seem to explain a lot
# alternative: selecting the number of components by having a threshold on pve, f.ex 80% or 70%

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

# Still not very clear, but the second sub-plot seems to show a small tendency towards 
# two groups..?

# Extra: hierarchical clustering *after* a PCA

# Hierarchical clustering on the SCORES (weights) of the first 3 PCs
data.dist.pca3 <- dist(pr.out$x[,1:3], method = 'manhattan')
hclust.df <- hclust(data.dist.pca3, method="complete" )
par(mfrow=c(1,1))
plot(hclust.df, labels=groups)
predicted <- cutree(hclust.df, k=2)
table(predicted, groups)

# K-means clustering on the SCORES (weights) of the first 3 PCs
predicted.kmean <- kmeans(pr.out$x[,1:3], 2, nstart=20)$cluster
table(predicted.kmean, groups)

## Possible comment / conclusion:
## -> this last clustering on the PCA is explaining 
## the Leukemia groups much better!

## Remember: you do NOT need to do all of these analyses on the exam! One is enough. 

## Final note: why are the methods not working so well on this dataset?? These types of results are 
## actually quite common, since real datasets are often quite messy and does not necessarily have
## a clear grouping structures. So: report what you see, even if the results are not "pretty".

