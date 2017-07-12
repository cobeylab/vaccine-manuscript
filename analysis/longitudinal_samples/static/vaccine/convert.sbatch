#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --partition=amd


rm results/$SLURM_ARRAY_TASK_ID/long.csv

if [ -f results/$SLURM_ARRAY_TASK_ID/out.long ]; then
    Rscript convert.R $SLURM_ARRAY_TASK_ID
fi

