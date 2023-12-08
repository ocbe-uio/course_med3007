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


# load the data and compute some summary statistics.
data.birds <- read.table('./data/birds.txt')
head(data.birds)
is(data.birds)

# the data look fine, and the data object is a data frame
apply(data.birds, 2, mean)
apply(data.birds, 2, sd)
apply(data.birds, 2, min)
apply(data.birds, 2, max)

summary(data.birds)

# plot the data in the way you think is best.
plot(data.birds, pch=16)

# use hierarchical clustering to prove the biologist's theory.
h.birds <- hclust(dist(data.birds, method='manhattan'), method='single')
plot(h.birds)
clust.birds <- cutree(h.birds, k=2)
# the two groups are pretty clear.. are they the correct ones?

plot(data.birds, pch=16, col = clust.birds)
# :)


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

# load in the data using read.csv(). You will need to select header=F.
data <- read.csv("./data/Ch10Ex11.csv", header=FALSE)
data <- t(data) # want each row to represent a sample ... should have n=40 samples/rows

# do hierarchical clustering and k-means

# hierarchical clustering
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

# kmeans clustering
predicted.kmean <- kmeans(data, 2, nstart=20)$cluster
table(predicted.kmean, true.groups )
# also very well!
