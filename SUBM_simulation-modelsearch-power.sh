#!/bin/bash

#$ -N simulation-modelsearch-power  # Job name
#$ -t 1:100    # Number of jobs
#$ -q long.q     # Queue. Use long.q for run time >8h all.q
#$ -l h_vmem=1G # Reserve 2 GB memory 
#$ -tc 128      # Max concurrent jobs
#$ -cwd         # Run in current directory
#$ -o output/simulation-modelsearch-power/   # Direct output to subdirectory
#$ -e output/simulation-modelsearch-power/   # Direct output to subdirectory

R CMD BATCH BATCH_simulation-modelsearch-power.R output/simulation-modelsearch-power/$JOB_NAME-I-$SGE_TASK_ID.Rout

##### cd Cluster/LVMproject/modelsearch
##### rm -r ./output/simulation-modelsearch-power/*
##### rm -r ./Results/simulation-modelsearch-power/*
##### qsub SUBM_simulation-modelsearch-power.sh

##### Your job-array 12746.1-100:1 ("simulation-modelsearch-power") has been submitted
##### 01/08/20 2:41 

##### qstat -u hpl802
## for n =  30, 300, 500
#     user   system  elapsed 
# 61444.58    37.35 


