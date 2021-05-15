#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=stdout.txt
#SBATCH --error=stderr.txt
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --partition=cobey

./vaccine_database.py
./write_params_to_sqlite.py

#rm -rf results
