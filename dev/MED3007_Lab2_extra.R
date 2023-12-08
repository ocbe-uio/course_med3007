###-----------------------------------------------------------------------------------###
###     Data visualization and dimensional reduction with Rstudio: extra material     ###
###-----------------------------------------------------------------------------------###


###------------------------###
### 1. MATRIX OBJECTS IN R ###
###------------------------###

### scalars
###--------

a <- 1
a

### vectors
###--------

# vectors as collections of numbers (quite stupid)
v <- c(2,3,5,4)
v
# vectors as sequences
u <- seq(2,5,len=3)
u
z1 <- seq(0,1,by=.01)
z1
z2 <- 1:5
z2

# vectors as repetitions
z <- rep(2, 4)
z
z1 <- rep(1:2, 4)
z1
z2 <- rep(1:2, each=4)
z2

### matrices
###---------

# create a matrix from a vector of numbers
W <- matrix(data = c(11,12,13,14,15,16), nrow = 2, ncol = 3, byrow = FALSE)
W
# create a matrix by merging different vectors (by row)
W <- rbind(c(11,13,15),c(12,14,16))
W
W <- cbind(c(11,13,15),c(12,14,16))
W
# (by column)
W <- cbind(c(11,12),c(13,14),c(15,16))
W


### take elements from vectors and matrices: [  ]
# vectors are indexed by a single number
v
v[2]
v[c(2,3)]
v[2:3]

# matrices are indexed by two numbers,
# the first for rows and the second for columns
W
W[2,3]
W[2,c(2,3)]
W[2,]
W[,c(2,3)]

### arrays
###-------

## arrays are just the generalization of a matrix to more that 2 dimensions..
## not really so used in practice, unless you need advanced statistics

myarray <- array(1:40, dim = c(4,5,2))
myarray
# for accessing elements: the same you do with matrices
myarray[,,1]
myarray[2,3,]

###------------------------------###
### 2. MATRIX & ARRAY OPERATIONS ###
###------------------------------###

## the easiest way to do matrix & array operations is with the function
## "apply" (and all its variations)

## you can use the function "apply" to apply a function on columns or rows
## of a matrix without using an explicit loop. Syntax is as follows:
## apply(X, margin, function)
## "X" is the array (matrix),
## "margin" is the dimension that should be kept
## "function" is the function that is to be applied for each component

## remark: margin needs to be an available dimension!

# example from the slides of the first day
mat <- matrix(c(1:6), ncol=2, byrow=TRUE)
mat 
apply(mat, 1, sum)
apply(mat, 2, sum)

# let's try with the array
myarray
apply(myarray, 3, sum)
apply(myarray, 2, sum)

## same commands exist for list objects: sapply/lapply
a1 <- 1:10                
a2 <- 2:30
a3 <- 3:40
lst <- list(a1,a2,a3)
lst

lapply(lst, mean)
sapply(lst, mean)




###---------------------------------------###
### SOME MORE EXPLORATION & VISUALIZATION ###
###---------------------------------------###

## Consider the dataset "food.txt", concerning the different kinds of food
## consumption in a number of European states.
## The purpose is to use this toy example to get used to PCA interpretation

food <- read.table('./data/Food.txt', header=T)
head(food)
n <- dim(food)[1]
p <- dim(food)[2]
n
p

# let us have a look at the variables

boxplot(food, main = 'Food consumption across European countries')
# quite different variability of the variables -> I scale the data in the PCA

## PCA
##----
pc.food <- prcomp(food, scale=TRUE)
pc.food
summary(pc.food)

## let us explore the PCA results

## proportion of explained variance (same as above!)
pc.food.var <- pc.food$sdev^2
pc.food.pve <- pc.food.var/sum(pc.food.var)
pc.food.pve
cumsum(pc.food.pve)

## loadings
# (coefficients in the linear combination of the variables
#  that defines each component)
load.food    <- pc.food$rotation
load.food


## scores
# (projections of each datum in the new space defined by the components)
# (if we select some of the components, dimensional reduction of the data)
scores.food <- pc.food$x
scores.food

## SIDE NOTE:
# loadings and PC variances are respectively
# eigenvectors and eigenvalues of the data correlation matrix
pc.food$rotation
pc.food$rotation[1:p,1:p]
pc.food.var
eigen(cor(food))


## Some PCA visualisation
##-----------------------

## 1. graphical summary of the proportion of variance explained
layout(matrix(c(2,3,1,3),2,byrow=T))
# variance of the PCs
barplot(pc.food.var, las=2, main='Principal components', ylab='Variances', 
        names.arg = paste('PC ',1:length(pc.food.var),sep=''))
# variance of the original data variables
barplot(apply(food, 2, sd)^2, las=2, main='Original Variables', ylim=c(0,150), ylab='Variances')
# PVE
plot(cumsum(pc.food.pve), type='b', axes=F, xlab='number of components',
     ylab='contribution to total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(food),labels=1:ncol(food),las=2)


## 2. dispersion on original data and in scores
layout(matrix(c(1,2),2))
boxplot(food, las=2, col='red', main='Original Variables')
scores.food <- data.frame(scores.food)
boxplot(scores.food, las=2, col='red', main='Principal components')

## 3. correlation among scores compared to original variables
pairs(food, asp=1)
pairs(scores.food, asp=1)


## 4. loadings of the first 3 PCs (those that might be the relevant ones)
par(mar = c(1,4,0,2), mfrow = c(3,1))
for(i in 1:3) barplot(load.food[,i], ylim = c(-1, 1))

## 5. a possible comparison: scatterplot of the data corresponding to the 2 variables
##                           with larger variance, and of the scores of the first
##                           2 principal components
layout(matrix(c(1,2),1))
plot(food[,'Cereals'],food[,'Milk'],type="n",xlab="Cereals",ylab="Milk", asp=1)
text(food[,'Cereals'],food[,'Milk'],dimnames(food)[[1]], cex=.7)
plot(scores.food[,1],-scores.food[,2],type="n",xlab="pc1",ylab="-pc2", asp=1)
text(scores.food[,1],-scores.food[,2],dimnames(food)[[1]], cex=.7)


