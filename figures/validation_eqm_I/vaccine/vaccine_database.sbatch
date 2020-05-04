#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=stdout.txt
#SBATCH --error=stderr.txt
#SBATCH --time=03:30:00
#SBATCH --partition=ivyb

./vaccine_database.py > database.txt
