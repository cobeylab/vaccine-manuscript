#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem=2000
#SBATCH --time=00:30:00
#SBATCH --partition=broadwl


let START=51
let END=51

# if end>total tasks then end

for SUB_JOB_NUM in `seq $START $END`
do
	rm results/$SUB_JOB_NUM/long.csv
	Rscript convert.R $SUB_JOB_NUM
done

# Use this instead to restart runs that go extinct:
# ${ANTIGEN_ROOT}/run_py --with-restarts
