###---------------------------------###
###     Clustering with Rstudio     ###
###---------------------------------###

### OUTLINE
###---------

### 1 Clustering: simple example
### 2. Exercise
### BREAK
### 3. Clustering: NCI60 Data
### 4. Exercise

## These R exercises are partly adopted from "James et al. (2013).
## An Introduction to Statistical Learning with Applications in R"
## (Chapter 10: "Unsupervised Learning").

# setwd("~/Dropbox_UiO/Dropbox/MED3007_2023/day 4")

###-------------------------------###
### 1. CLUSTERING: SIMPLE EXAMPLE ###
###-------------------------------###

## Consider the IRIS dataset, where we *have* the data labels, but we forget about them, and use them
## for validating / interpreting the clusters.

## We first explore the data, and see how the metric choice
## affects the distance matrix

species.name <- iris[,5]
iris4 <- iris[,1:4]
n <- dim(iris4)[1]
p <- dim(iris4)[2]

# let us plot the data and remember about them
pairs(iris4)

# let's do some jittering (usually done for this dataset)
iris4 <- iris4 + cbind(rnorm(n, sd=0.025))
pairs(iris4)

# we now want to compute the distance matrix of the iris data
# we choose the Euclidean metric first, and then we vary this to see the effect
iris.e <- dist(iris4, method='euclidean')
is(iris.e)
help(dist)

par(mfrow=c(1,1))
image(as.matrix(iris.e), main='Euclidean metric', xlab='i', ylab='j')
# REMARK: this is the distance matrix resulting from a "toy" example!!!
# in reality, data are obviously never ordered according to true labels.
# So, let's re-arrange the data to have a more realistic image.

myshuffle <- sample(n)
iris4 <- iris4[myshuffle,]
iris.e <- dist(iris4, method='euclidean')

# let us try and see what happens when using a different metric
iris.m <- dist(iris4, method='manhattan')
iris.c <- dist(iris4, method='canberra')

par(mfrow=c(1,3))
image(as.matrix(iris.e), main='Euclidean metric', xlab='i', ylab='j')
image(as.matrix(iris.c), main='Canberra metric', xlab='i', ylab='j')
image(as.matrix(iris.m), main='Manhattan metric', xlab='i', ylab='j')



## We now try and perform the hierarchical clustering of the iris dataset

# with Euclidean metric
iris.es <- hclust(iris.e, method='single')
iris.ea <- hclust(iris.e, method='average')
iris.ec <- hclust(iris.e, method='complete')

# with  Canberra metric
iris.cs <- hclust(iris.c, method='single')
iris.ca <- hclust(iris.c, method='average')
iris.cc <- hclust(iris.c, method='complete')

# Remark: if you'd like to have detailed information for ex on hierarchical
#         clustering with euclidean-complete combination, just have a look
#         at the corresponding R object (iris.ec)
# --> More info in the script "MED3007_Lab3_extra.R".

# plot all dendrograms
# (Note how "smart" R is: the plot command automatically recognizes the object
#  of class "hclust", and plots the dendrogram without needing to specify)
par(mfrow=c(2,3))
plot(iris.es, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(iris.ec, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(iris.ea, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(iris.cs, main='canberra-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(iris.cc, main='canberra-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
plot(iris.ca, main='canberra-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')


# The choice of linkage certainly affects the results, and it does much more than
# the metric choice: even though we also see some differences between euclidean and
# canberra results, the metric does not change the shape of the dendrogram in the
# same way as the linkage does.

# Typically, single linkage will tend to yield trailing clusters:
# very large clusters onto which individual observations attach one-by-one.

# On the other hand, complete linkage is often preferred
# because it results in the most compact clusters.

## Let us now add some rectangles for actually dividing the data into clusters,
## and in particular 3 clusters

# we can do this by using the function "rect.hclust"
par(mfrow=c(2,3))
plot(iris.es, main='euclidean-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(iris.es, k=2)
rect.hclust(iris.es, k=3)
plot(iris.ec, main='euclidean-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(iris.ec, k=2)
rect.hclust(iris.ec, k=3)
plot(iris.ea, main='euclidean-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(iris.ea, k=2)
rect.hclust(iris.ea, k=3)
plot(iris.cs, main='canberra-single', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(iris.cs, k=2)
rect.hclust(iris.cs, k=3)
plot(iris.cc, main='canberra-complete', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(iris.cc, k=2)
rect.hclust(iris.cc, k=3)
plot(iris.ca, main='canberra-average', hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
rect.hclust(iris.ca, k=2)
rect.hclust(iris.ca, k=3)

# Complete linkage is clearly the best.
# The dendrograms also quite clearly show the presence of the 3 groups.
# The Euclidean distance with complete linkage seems to give the best result,
# but the Canberra distance could have worked as well.

## In order to decide between the Euclidean and Canberra distance, let us
## compare how coherent the 3 groups are in the two cases.

# we use the function "cutree" to "cut" the dendrogram at a given number of groups
cluster.ec <- cutree(iris.ec, k=3)
cluster.cc <- cutree(iris.cc, k=3)
# we have to compare the two clustering results with each other.
# we can use the function "table"
table(cluster.cc, cluster.ec)
# both cases completely agree on one cluster (this we see from the plots also),
# while they disagree somewhat more in the other two clusters.

###--------------###
### 2. EXERCISE  ###
###--------------###

## In the dataset "birds.txt" are included the measurements in cm of the
## length and width of the chest of 50 birds of a specific species (Passer domesticus).
## The biologist who has taken the measurements wants to prove that this kind of
## birds can be grouped  in 2 clusters from their chest characteristics.
## Help him to prove his theory.
## 1. Load the data and compute some summary statistics.
## 2. Plot the data in the way you think is best.
## 3. Use hierarchical clustering to prove the biologist's theory.



###------------------------------###
### 3. CLUSTERING THE NCI60 DATA ###
###------------------------------###

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

## Clustering the Observations of the NCI60 Data via K-means clustering
##---------------------------------------------------------------------

## Run K-means clustering with K=2,3, and 4 clusters.
## Why do we have to set a random seed (set.seed())?

# K-Means Clustering with K=2, 3, 4 clusters:
set.seed(4)
km.out2 <- kmeans(sd.data, 2, nstart=20)
km.out3 <- kmeans(sd.data, 3, nstart=20)
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


###-------------###
### 4. EXERCISE ###
###-------------###


## Consider again the gene expression data set "Ch10Ex11.csv"
## (which can be also found on the book website, www.StatLearning.com)
## that consists of 40 tissue samples with measurements on 1,000 genes.
## The first 20 samples are from healthy patients,
## while the second 20 are from a diseased group.

## 1. Load in the data using read.csv(). You will need to select header=F.
##    Alternatively: load in the data using "Import dataset" in the upper
##    right window and click "no" in the "Header" option. You also need to
##    transform the dataset using t().
## 2. Perform hierarchical clustering and k-means on these data,
#     and compare the results. Can the original groups be found?
