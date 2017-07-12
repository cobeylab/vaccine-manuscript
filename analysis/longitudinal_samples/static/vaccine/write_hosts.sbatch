#!/bin/bash

#SBATCH --ntasks-per-node=1
#SBATCH --output=stdout.txt
#SBATCH --error=stderr.txt
#SBATCH --time=00:30:00
#SBATCH --partition=amd
#SBATCH --mem-per-cpu=3000

Rscript write_hosts.R 1 static_10percent.csv
box upload static_10percent.csv vaccine_evolution/static_10percent.csv
