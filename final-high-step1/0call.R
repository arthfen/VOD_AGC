################### Script: 1call.R ##################################
## This script starts the procesing cluster and define the basis dimensions
## It follows by calling the next scripts, '1smoothers.R' and '2fit_ML.R'

## Creates the parallel processing cluster of 'nodes' nodes
library(snow)
nodes <- 12
cl <- makeCluster(nodes, type = 'MPI')

## Sets a new environment for the scripts to be run
nenv = new.env()

## Creates the dataframe with the basis dimensions for each smoother
### In order to understand this, go to file 2fit_ML.R
# From b5 on, values are lower due to the small n. of unique values for TMF shares and for the 'year' variable
df <- data.frame(b0 = 12, b1 = 12, b2 = 12, b3 = 12, b4 = 12, b5 = 4, b6 = 8)
df$res <- 0

h.b0 = df[1, 'b0']
h.b1 = df[1, 'b1']
h.b2 = df[1, 'b2']
h.b3 = df[1, 'b3']
h.b4 = df[1, 'b4']
h.b5 = df[1, 'b5']
h.b6 = df[1, 'b6']
h.sp <- 'cs' # We are using cubic splines (see the documentation of the mgcv package)
h.m <- 3

source('1smoothers.R', local = nenv)
source('2fit_ML.R', local = nenv)
