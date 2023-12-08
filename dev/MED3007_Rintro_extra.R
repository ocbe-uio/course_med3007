###--------------------------------------------------###
###     Introduction to RStudio: extra material      ###
###--------------------------------------------------###


###------------------------###
### 1.1 BASIC OBJECTS IN R ###
###------------------------###
### Most of these objects are needed to handle the data/results in R

### lists:
### these objects are VERY IMPORTANT because typically results of analyses
### performed with R functions have a list structure
### as the name says, they are list of objects of any kind

# they can be created with the function "list"
course <- list(title  = 'Statistical Genomics',  
               date  = '12/01/2021',
               numb_hours = 25,
               outcome_data = 10:1,
               covariate_data =  matrix(1:40,10,4)
)
course

# or they can be returned by an R function performing something (example: t test)
data <- rnorm(100,1) # 100 random numbers generated from normal distribution with mean=1 and sd=1
result <- t.test(data,alternative="two.sided")
names(result)     # check the list content (only if the elements have names)
result$p.value
result$conf.int

# some useful functions
length(course)
length(result)
# element extraction
result[[3]]
course[[2]]
course$date
M1 <- course$covariate_data # M1 is a matrix
dim(M1)

### dataframes: these are "special types" of list 
###             all elements MUST be vector of the same length
###             this is the form of data when they are IMPORTED into R (from txt, csv, ..)
###             (see later, when "importing data into R")

# they can be created with the function data.frame, even though this is less frequent than
# getting a dataframe from imported data (same as for the function "list")
data_matrix <- data.frame(cbind(course$outcome_data, course$covariate_data))
names(data_matrix) <- c('outcome','V1','V2','V3','V4')
data_matrix
# more frequently
data_matrix <- read.table('./data/CHDAGEdata.txt',header=TRUE)
head(data_matrix) # useful function for inspecting the top of a dataframe
is(data_matrix) # which class does the object come from? 
# hint: you can also convert loaded data into a dataframe by using the function as.data.frame


# CAUTION: dataframes look like matrices, but they are not!
#          however, some matrix commands can be used (this will be clear tomorrow)
data_matrix[2,3]
dim(data_matrix)
# I add a new column, for example
gender <- c(rep('M',50),rep('F',50))
data_matrix <- cbind(data_matrix, GENDER = as.factor(gender))
head(data_matrix)

# IMPORTANT: a dataframe is an environment. In R this means that a dataframe is closed
# with respect to the global environment; thus, the variables INSIDE the dataframe
# cannot be accessed from the global environment
AGE
# but, they can be accessed if the dataframe becomes "attached"
attach(data_matrix)
AGE
# CAUTION: this can mask other variables which live in the global environment
#          (if any variable inside the dataframe has the same name as a variable in the global environment)
detach(data_matrix)

# when the dataframe is not attached, the variables inside it
# can always be accessed as it is done for lists
data_matrix$AGE
data_matrix$GENDER
# OR as it is done for matrices
data_matrix[,2]
data_matrix[,4]
# OR with variables names
data_matrix[,'AGE']
data_matrix[,'GENDER']


### PACKAGES
###---------

### Packages are collections of R functions, data, and compiled code,
### that can be loaded and directly used in R.
### The material is collected in a well-defined format, so that for each
### element there is a help page to guide its use.

# The directory where packages are stored is called the library.
# To load a package in R it is thus sufficient to use the command "library"

library(MASS)
# "MASS" is a package that contains many interesting data
data()
# if you want to load a particular dataset you just add it as an argument to "data"
data(AirPassengers)

# many packages (like MASS) are already installed in R
# but most packages need to be installed first, and then they can be loaded when needed

# ex: package for generalized additive models
install.packages('gam') # this is needed only the first time
library(gam)            # this is needed every time you open R and need this package


### HELP
###-----

# if you know the function name
# help(FUNCTIONNAME) OR ?FUNCTIONNAME
help(hist)
?plot

# if you don't remember the function name
help.start()

# OR in RStudio: "Help" --> "R Help")
# OR google! :)



