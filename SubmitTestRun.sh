#!/bin/bash

#SBATCH -t 1- 
#SBATCH -N 1 
#SBATCH -n 1

# Ensure Rlogs directory exists
mkdir -p ./Rlogs

# Load the R module
module add r/4.1.3

R CMD BATCH --vanilla --args --shortname=Test --core=1 --copy=1 --GibbsRun=T --algorithm=DP --PostD=4  --L=1 --folder=Tests --n=100 --user=emmamit TestRun.R ./Rlogs/TestRun.out
