# code for clustering

# food ----
# food example for clustering
# the clustered groups are easier to interpret, hence we use it for teaching

# load data, and change the title of pulse to Legume
food <- read.table('./lab/data/Food.txt', header=T)
colnames(food)[7] <- 'Legume'
head(food) # first 6 lines of code

# the food data has columns (features) that have quite big difference
# Scale the data to zero mean and unit variance, call it food_s
food_s <- scale(food)


## hierarchical clustering ----

# HC requires distances (not original data) as input
# Calculate the distance matrix on the scaled data
# (default method for distance = Euclidean)
food_dist <- dist(food_s)
# food_dist <- dist(food_s, method="euclidean")
# food_dist

# food_dist is a vector of 300 elements
length(food_dist) # 300

# this is the number of pairs of data
# we have 25 data points, one for each country
# Albania vs Austria, Albania vs Belg.Lux, ...
# in total you have 25 * 24 / 2 pairs



# Perform clustering with different linkage methods:
hc.complete <- hclust(food_dist, method="complete")
plot(hc.complete, labels=rownames(food), main="Complete Linkage", xlab="", sub="")

# single linkage
hc.single <- hclust(food_dist, method="single")
plot(hc.single, labels=rownames(food), main="Single Linkage", xlab="", sub="")

# average linkage
hc.average <- hclust(food_dist, method="average")
plot(hc.average, labels=rownames(food), main="Average Linkage", xlab="", sub="")


# unscaled data, complete linkage
hc.unscaled <- hclust(dist(food), method="complete")
plot(hc.unscaled, labels=rownames(food), main="Complete linkage with unscaled features", xlab="", sub="")


# correlation on scaled data, complete linkage
cor(food_s) # this is by col (food)
cor(t(food_s)) # by country

# 1 minus makes them all positive
dd <- as.dist(1-cor(t(food_s)))
hc.corr <- hclust(dd, method="complete")
plot(hc.corr, labels=rownames(food), main="Complete linkage with correlation-based distance", xlab="", sub="")


## heatmap ----
# heatmap (default)
# the default heatmap does clustering for both row and col variables
# dendrograms are also shown
heatmap(food_s)

# you can specify the arguments so that no clustering is done
# the original order of col and row are kept
# heatmap with no clustering
heatmap(food_s, Rowv = NA, Colv = NA)

# can also do clustering for only row (or col)
heatmap(food_s, Colv = NA)




# _________ ----
# NCI 60 ----
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

# hierarchical clustering ----

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
par(mfrow=c(1,1))
plot(hc.complete, labels=nci.labs, main="Complete Linkage", xlab="", sub="")
# plot(hc.average, labels=nci.labs, main="Average Linkage", xlab="", sub="")
# plot(hc.single, labels=nci.labs,  main="Single Linkage", xlab="", sub="")

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







# k means ----
## Clustering the Observations of the NCI60 Data via K-means clustering

## Run K-means clustering with K=2,3, and 4 clusters.
## Why do we have to set a random seed (set.seed())?

# K-Means Clustering with K=2, 3, 4 clusters:
set.seed(4)
# km.out2 <- kmeans(sd.data, 2, nstart=20)
# km.out3 <- kmeans(sd.data, 3, nstart=20)
km.out4 <- kmeans(sd.data, 4, nstart=20)

## Read the help file ?kmeans to understand what the argument nstart=20 does.
## Comparing an analysis with nstart=20 versus nstart=1 demonstrates
## how the cluster results can be improved if we allow more evaluations
## with different randomly chosen starting centroids.

?kmeans

# More evaluations with different starting centroids improve the clustering:
set.seed(3)
km.out <- kmeans(sd.data, 3, nstart=1)
km.out$tot.withinss

km.out <- kmeans(sd.data, 3, nstart=20)
km.out$tot.withinss

## Next, let us compare the K-means and hclust solutions with 4 clusters.
## Interpret the results.

# first let's look at the clustering result with 4 groups
km.out4$cluster

# then, we can directly compare the k-means result (along rows)
# with the hierarchical clustering result (along columns)
table(km.out4$cluster, hc.clusters[,"4"], deparse.level=2)

# Quite a bit of disagreement between the two methods, except for maybe two clusters.


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



# heatmap ----
# find out how to set argument for heatmap by reading the documentation
?heatmap

# We use the scores of the PCA on the NCI60 data, to reduce dimension

# default heatmap (clustering done on both col and row)
heatmap(pr.out$x)

# clustering (on correlation) done for rows
heatmap(pr.out$x, Rowv = as.dendrogram(hc.corr), Colv = NA)

# I now plot less components for the sake of clarity,
# I add tumor type, and I give a title
par(cex.main = .7)
heatmap(pr.out$x[,1:40], Rowv = as.dendrogram(hc.corr), Colv = NA,
        labRow = nci.labs, main = 'Heatmap of the scores of the first 40 PCs on the NCI60 data')







# ________ ----
# CH10Ex11 ----

# load in the data using read.csv(). You will need to select header=F.
# set the right path to load the data!
data <- read.csv("lab/data/Ch12Ex13.csv", header=FALSE)
data <- t(data) # want each row to represent a sample ... should have n=40 samples/rows

# do hierarchical clustering and k-means

# hierarchical clustering ----
data.dist <- dist(data) # need to compute the distance matrix
hclust.df <- hclust(data.dist, method="complete" )
#alternatives:
#hclust.df <- hclust( D, method="average" )
#hclust.df <- hclust( D, method="single" )

# find the clusters
predicted <- cutree( hclust.df, k=2 )
true.groups <- c( rep(0,20), rep(1,20) )

# How well does our clustering predict health vs. diseased
table(predicted, true.groups )
# very well!!

# kmeans  ----

predicted.kmean <- kmeans(data, 2, nstart=20)$cluster
table(predicted.kmean, true.groups )
# also very well!






