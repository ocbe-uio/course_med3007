###---------------------------------###
###     Introduction to RStudio     ###
###---------------------------------###

### OUTLINE
###---------

### 1. RStudio: how does it work
### 2. Make a script (a file with code and comments)
### 3. Load in a dataset
### 4. Exercise


# Set the working directory:
# This can be done by going to Session -> Set working directory -> To source file location 
# .. or by writing a code
getwd() # get source file location path
setwd("/Users/emilie.odegaard/Dropbox_UiO/Dropbox/MED3007_2023/day 1")

# Why do we need to state this? 
# We often need to specify the folder where R is going to look for files (ex. data),
# or save files (ex. plots, results).


###------------------------###
### 1.2 DATA IMPORT/EXPORT ###
###------------------------###
### The most common way of importing data in R is loading them from a file

### DATA IMPORT
###------------

### from TXT or DAT files
###----------------------

mydata.pima <- read.table('./data/pima.txt', header=TRUE)
# exaclty the same with .DAT format!
# let's look at the data
mydata.pima
# too many lines! with this function, you only check the top of the data frame
head(mydata.pima)

### from a CSV file
###----------------

# if the first row contains the variable names, comma is the separator
# and I want to assign the variable id to row names:
mydata.listeria <- read.csv('./data/listeria.csv', header=TRUE, sep=",")
head(mydata.listeria)
dim(mydata.listeria)

### from EXCEL
###-----------

## Using the RStudio interface: in the top-right panel, in the "Environment" tab,
##    click "Import Dataset" and then "From Excel". The data will also appear in the "View" mode
## Remark: this is EQUIVALENT  to the following lines
install.packages('readxl')
library(readxl)
mydata.excel <- read_excel("data/Testfil_Rcourse.xlsx")
View(mydata.excel)
head(mydata.excel)
# Remark: mydata.excel is not technically a data.frame now, but a tibble!
#         This is a more difficult object to handle (even though the "tidyverse
#         approach" has advantages; read about the tidyverse R package)
#         Hence, better to transform to a data.frame
mydata.excel <- as.data.frame(mydata.excel)

### from R data objects
###--------------------
# This might happen when you have previously performed some analyses and would like to 
# load the results in R. Or, if someone else prepared the data for you by doing some
# pre-processing in R (typical with genomics data)

load('./data/mydata.RData')
# RData objects are collections of any R-friendly objects (as many as you like)
# that can be loaded directly into your working session


### DATA EXPORT ###
###-------------###

## how to save TO TXT
write.table(mydata, file = './data/mydata.txt')
# CAUTION: the file "mydata.txt" containing the test data is saved in the working directory
#          previously selected (if selected), otherwise in the default directory 
#          (NEVER use the default directory!!)

# to save different objects in the same file
save(mydata.excel, mydata.pima, mydata.listeria, file = './data/various_objects.RData')
# a .RData file is an R workspace, which can contain as many variables as needed
# if you double click it, R is directly opened
# (you can have many Rs opened at the same time..)


###----------------###
### 1.3 EXERCISE 1 ###
###----------------###

# Load the data file "testRcourse.txt" and access the variable "age"
# Draw a histogram of this variable and save it
# Remember to create a script!





