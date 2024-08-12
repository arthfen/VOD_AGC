################### Script: 1call.R ##################################
## This script starts the procesing cluster and follows by calling the next script, '2fit_ML.R'

## Creates the parallel processing cluster of 'nodes' nodes
library(parallel)
nodes <- 12
cl <- makeCluster(nodes, type = 'MPI')

## Sets a new environment for the scripts to be run
nenv = new.env()
source('2fit_ML.R', local = nenv)

