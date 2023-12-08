###-------------------------------------------------###
###     Clustering with Rstudio: extra material     ###
###-------------------------------------------------###

## Clustering the Observations of the NCI60 Data via Hierarchical clustering
##--------------------------------------------------------------------------

library(ISLR)
nci.labs <- NCI60$labs # Sample labels (tissue type)
nci.data <- NCI60$data # Gene expression data set

## We start by scaling the data, calculate the distance matrix
## (using the Euclidean distance), and then investigate different linkage methods.

# Scale the data to zero mean and unit variance:
sd.data <- scale(nci.data)

# Calculate the distance matrix (default = Euclidean):
data.dist <- dist(sd.data)
data.dist <- dist(sd.data, method="euclidean")

# Perform clustering with different linkage methods:
hc.complete <- hclust(data.dist, method="complete")
hc.average <- hclust(data.dist, method="average")
hc.single <- hclust(data.dist, method="single")

names(hc.complete)
hc.complete$merge  # order of aggregations of samples / clusters
hc.complete$height # distance at which aggregations happen
hc.complete$order  # correct order of the samples for obtaining the plot
hc.complete$labels # labels (numeric, since we don't know the original categories!)
hc.complete$method
hc.complete$call
hc.complete$dist.method

help(hclust)
is(hc.complete)

## Compare the dendrograms and comment on the results in light of how
## the different linkage methods work.
par(mfrow=c(1,3))
plot(hc.complete, labels=nci.labs, main="Complete Linkage", xlab="", sub="")
plot(hc.average, labels=nci.labs, main="Average Linkage", xlab="", sub="")
plot(hc.single, labels=nci.labs,  main="Single Linkage", xlab="", sub="")

# In light of what we see, we continue by further investigating
# the complete linkage clustering.

## First, we use cutree() to compare the results when the data are separated
## into either 2 or 4 clusters.

# Compare 2 clusters and 4 clusters:
hc.clusters <- cutree(hc.complete, c(2, 4))
table(hc.clusters[,"2"], hc.clusters[,"4"])

# How are the labels distributed between clusters:
table(hc.clusters[,"4"], nci.labs)

# visualize the cuts
par(mfrow=c(1,1))
plot(hc.complete, labels=nci.labs, main="Complete Linkage", xlab="", sub="")
abline(h=140, col="red")  # 4 clusters
abline(h=150, col="blue") # 2 clusters

## Finally, we see what happens if we use unscaled data instead of scaled data,
## or if we use a correlation-based distance metric instead of the Euclidean distance.
## Compare the dendrograms:
## How different are the resulting clusterings?
## Do you recognise subclusters that are consistent?


# Compare scaled data versus non-scaled data:
hc.unscaled <- hclust(dist(nci.data), method="complete")
par(mfrow=c(1,1))
plot(hc.unscaled, labels=nci.labs, main="Complete linkage with unscaled features", xlab="", sub="")

# Compare Euclidean distance with correlation-based distance:
dd <- as.dist(1-cor(t(sd.data)))
hc.corr <- hclust(dd, method="complete")
par(mfrow=c(1,1))
plot(hc.corr, labels=nci.labs, main="Complete linkage with correlation-based distance", xlab="", sub="")


## Clustering the Observations of the NCI60 Data via K-means clustering
##---------------------------------------------------------------------

## Run K-means clustering with K=2,3, and 4 clusters.
## Why do we have to set a random seed (set.seed())?

# K-Means Clustering with K=2, 3, 4 clusters:
set.seed(4)
km.out2 <- kmeans(sd.data, 2, nstart=20)
km.out3 <- kmeans(sd.data, 3, nstart=20)
km.out4 <- kmeans(sd.data, 4, nstart=20)


## We can visualise the K-means clustering results of high-dimensional
## data by using PCA for dimension reduction first. We plot the first
## two principal components and colour the data points
## (= individual cell lines) by their assigned cluster.

# first, run PCA again on the NCI60 data
pr.out <- prcomp(nci.data, scale=TRUE)

# we can now visualise the K-Means results by labelling the data points
# in a plot of the scores of the first 2 principal components:
par(mfrow=c(1,3))
plot(pr.out$x[,1:2], col=(km.out2$cluster+1), main="K-Means with K=2",
     xlab="PC 1", ylab="PC 2", pch=20)
plot(pr.out$x[,1:2], col=(km.out3$cluster+1), main="K-Means with K=3",
     xlab="PC 1",  ylab="PC 2", pch=20)
plot(pr.out$x[,1:2], col=(km.out4$cluster+1), main="K-Means with K=4",
     xlab="PC 1", ylab="PC 2", pch=20)

## heatmap
##--------

?heatmap

## let create a small heatmap, to fix ideas.
## We use the scores of the PCA on the NCI60 data, to reduce dimension

#  default choices
heatmap(pr.out$x) 

# I use the previous dendrogram for better ordering of the patients,
# and I remove the dendrogram for the components
heatmap(pr.out$x, Rowv = as.dendrogram(hc.corr), Colv = NA)

# I now plot less components for the sake of clarity,
# I add the patient's tumor type, and I give a title
par(cex.main = .7)
heatmap(pr.out$x[,1:40], Rowv = as.dendrogram(hc.corr), Colv = NA,
        labRow = nci.labs, main = 'Heatmap of the scores of the first 40 PCs on the NCI60 data')
# very nice plot!


## elbow plot
##-----------

names(km.out2)
# the within cluster sum-of-squares is within the object "withinss"
# but we need to run much more k-means in order to decide...

which.clust <- 1:15
within.clust.var <- NULL
for(k in which.clust){
  myresult <- mean(kmeans(sd.data, k, nstart=10)$withinss)
  within.clust.var <- c(within.clust.var, myresult)
}
# let's plot the values and look for the elbow
par(mfrow=c(1,1))
plot(which.clust, within.clust.var, type = 'b', lwd = 2, 
     xlab = 'number of clusters',
     ylab = 'within-cluster sum-of-squares', 
     main = 'k-means clustering of NCI60 data')

