#!/bin/bash

#$ -N simulation-modelsearch-type1  # Job name
#$ -t 1:100    # Number of jobs
#$ -q long.q     # Queue. Use long.q for run time >8h all.q
#$ -l h_vmem=6G # Reserve 2 GB memory 
#$ -tc 128      # Max concurrent jobs
#$ -cwd         # Run in current directory
#$ -o output/simulation-modelsearch-type1/   # Direct output to subdirectory
#$ -e output/simulation-modelsearch-type1/   # Direct output to subdirectory

R CMD BATCH BATCH_simulation-modelsearch-type1.R output/simulation-modelsearch-type1/$JOB_NAME-I-$SGE_TASK_ID.Rout

##### cd Cluster/LVMproject/article-multipleComparisons/
##### rm -r ./output/simulation-modelsearch-type1/*
##### rm -r ./Results/simulation-modelsearch-type1/*
##### qsub SUBM_simulation-modelsearch-type1.sh


##### Your job-array 12740.1-100:1 ("simulation-modelsearch-type1") has been submitted
##### 01/02/20 4:35 

##### qstat -u hpl802
#####      user    system   elapsed 
##### 10282.140    10.068 10286.672 
