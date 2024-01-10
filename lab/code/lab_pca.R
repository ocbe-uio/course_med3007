# code for EDA and PCA

# food  -----

food <- read.table('./lab/data/Food.txt', header=T)
head(food)
# pulses: legume
colnames(food)[7] <- 'Legume'

## explore and viz ----
# let us have a look at the variables
dim(food)

colnames(food)
rownames(food)

# subsetting
food[rownames(food) == 'Norway',]
food[rownames(food) %in% c('Norway', 'Denmark', 'Sweden'),]

# explore one variable
mean(food$Fish)
max(food$Fish)
min(food$Fish)

## visualize ----
# par(mfrow = c(1, 1))
# reset margin
# par(mar = c(5.1, 4.1, 4.1, 2.1))

hist(food$Fish)
hist(food$Fish, breaks = 20)

hist(food$Meat)


# boxplot
boxplot(food, las = 2,
        main = 'Food consumption across European countries')


# two variables
plot(food$Fish, food$Meat)

# add country
text(food$Fish, food$Meat, labels = rownames(food))

# choose aother pair of food
plot(food$Cereals, food$Fruit)
text(food$Cereals, food$Fruit, labels = rownames(food))

# pair-wise scatter
plot(food)


## PCA ----

pc_food <- prcomp(food, scale=TRUE)
# pc_food
summary(pc_food)



## loadings
# (coefficients in the linear combination of the variables
# that defines each component)
loading_food <- pc_food$rotation
loading_food
round(loading_food, digits = 2)

## scores
# (projections of each datum in the new space defined by the components)
# (if we select some of the components, dimensional reduction of the data)
scores_food <- pc_food$x
scores_food



## proportion of explained variance
pc_food_var <- pc_food$sdev^2
pc_food_pve <- pc_food_var/sum(pc_food_var)
pc_food_pve
cumsum(pc_food_pve)


## Some PCA visualisation
# par(mfrow = c(1,1))

# biplot
biplot(pc_food)
# biplot(pc_food, choices = c(1,2)) # first two pc

biplot(pc_food, choices = c(1,3)) # pc1, pc3



# graphical summary of the proportion of variance explained

# variance of the PCs
barplot(pc_food_var, las=2,
        main='Principal components', ylab='Variances',
        names.arg = paste('PC ',1:length(pc_food_var),sep=''))

# variance of the original data variables
barplot(apply(food, 2, sd)^2, las=2, main='Original Variables', ylim=c(0,150), ylab='Variances')


# PVE
plot(cumsum(pc_food_pve), type='b', axes=F, xlab='number of components',
     ylab='contribution to total variance', ylim=c(0,1))
abline(h=1, col='blue')
abline(h=0.8, lty=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(food),labels=1:ncol(food),las=2)




## OPTIONAL
# loadings and PC variances are respectively
# eigenvectors and eigenvalues of the data correlation matrix
pc.food$rotation
pc.food$rotation[1:p,1:p]
pc.food.var
eigen(cor(food))



# NCI60  ----

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





# Ch12Ex13 ----

# set the path to your own!
exp.data <-  read.csv("./lab/data/Ch12Ex13.csv",header=FALSE)

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


