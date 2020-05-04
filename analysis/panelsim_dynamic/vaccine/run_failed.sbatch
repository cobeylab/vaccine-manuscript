#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=29500
#SBATCH --time=03:00:00
#SBATCH --partition=broadwl


#echo $JOB
JOB=$SLURM_ARRAY_TASK_ID
cd $JOB
rm out.*
#rm parameters.json
rm parameters_out.json
rm vaccine.timeseries
rm samples.sqlite
rm hosts.tsv
rm output.sqlite
rm phenotypes.tsv
#mv parameters_out.json parameters.json

${ANTIGEN_ROOT}/run.py > "out.bashscreen"
excessDiversityFile="out.tmrcaLimit"
extinctFile="out.extinct"
../../process_job.py $JOB
cd ../

# Use this instead to restart runs that go extinct:
# ${ANTIGEN_ROOT}/run_py --with-restarts
