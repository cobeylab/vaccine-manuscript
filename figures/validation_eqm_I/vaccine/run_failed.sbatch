#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mem=24000
#SBATCH --time=1:00:00
#SBATCH --partition=sandyb

let SUB_JOB_NUM=$SLURM_ARRAY_TASK_ID

# if end>total tasks then end
cd $SUB_JOB_NUM
rm out.*
#rm parameters.json
rm vaccine.timeseries
rm samples.sqlite
rm hosts.tsv
rm output.sqlite
rm phenotypes.tsv

#mv parameters_out.json parameters.json
${ANTIGEN_ROOT}/run.py > "out.bashscreen"
excessDiversityFile="out.tmrcaLimit"
extinctFile="out.extinct"
../../process_job.py $SUB_JOB_NUM
cd ../

# Use this instead to restart runs that go extinct:
# ${ANTIGEN_ROOT}/run_py --with-restarts

