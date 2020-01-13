#!/bin/bash

#$ -N simulation-selection-type1  # Job name
#$ -t 1:50     # Number of jobs
#$ -q long.q    # Queue. Use long.q for run time >8h and all.q otherwise
#$ -l h_vmem=8G # Memory limit, e.g. reserve 1 GB memory 
#$ -tc 128      # Max concurrent jobs
#$ -cwd         # Run in current directory
#$ -o output/simulation-selection-type1/   # Direct output to subdirectory
#$ -e output/simulation-selection-type1/   # Direct output to subdirectory

R CMD BATCH BATCH_simulation-selection-type1.R output/simulation-selection-type1/$JOB_NAME-I-$SGE_TASK_ID.Rout --no-restore --no-save

## go to directory    ## cd Cluster/LVMproject/article-multipleComparisons/
## clean outputs      ## rm -r ./output/simulation-selection-type1/*
## clean results      ## rm -r ./Results/simulation-selection-type1/*
## submission command ## qsub SUBM_simulation-selection-type1.sh

## submission output  ## Your job-array 12745.1-50:1 ("simulation-selection-type1") has been submitted
## submission time    ## 01/08/20 2:31 

## documentation      ## https://ifsv.sund.ku.dk/biostat/biostatwiki/index.php/IT:Cluster
                      ## http://bayes/ganglia : current load and history can be monitored graphically at
                      ## qstat -f -u \*       : view the cluster queue type (the first column shows the job id)
                      ## qstat -j 1034        : show details of a job (or job array) with job id 1034 type     
                      ## qdel 1034            : delete the job with job id 1034 from the queue type
