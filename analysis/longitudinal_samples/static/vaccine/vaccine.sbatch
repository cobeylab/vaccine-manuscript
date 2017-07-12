#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem=28000
#SBATCH --time=04:00:00
#SBATCH --partition=amd

let START=$SLURM_ARRAY_TASK_ID*$N_PER_JOB
let END=$START+${N_PER_JOB}-1

# if end>total tasks then end

for SUB_JOB_NUM in `seq $START $END`
do
	cd $SUB_JOB_NUM
	${ANTIGEN_ROOT}/run.py > "out.bashscreen"
	excessDiversityFile="out.tmrcaLimit"
	extinctFile="out.extinct"
	../../process_job.py $SUB_JOB_NUM
	../../Rscript convert.R $SUB_JOB_NUM
	cd ../
done

# Use this instead to restart runs that go extinct:
# ${ANTIGEN_ROOT}/run_py --with-restarts
