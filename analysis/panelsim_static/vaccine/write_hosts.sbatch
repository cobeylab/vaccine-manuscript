#!/bin/bash

#SBATCH --ntasks-per-node=1
#SBATCH --output=stdout.txt
#SBATCH --error=stderr.txt
#SBATCH --time=01:30:00
#SBATCH --partition=cobey
#SBATCH --mem-per-cpu=8000

Rscript write_hosts.R 1 static_10percent_breadth_1_022820.csv
