#!/usr/bin/env python

import csv
import os
import json
import numpy as np
import subprocess
from collections import OrderedDict

DRY = False
#JOBS = [75,71,3,99,78,73,76,86,94,80,18,87,88,91,69,62,2,67,64,1,70,0,77,4,65,5,85,96,82,61,83,66,90,72,79,95,84,92,97,74]
JOBS = [51]
N_PER_JOB=1

def submitJob(rootDir, sweepDir, resultsDir):
	# Construct SLURM command
	sbatchFilename = os.path.join(sweepDir, 'run_failed.sbatch')

 	submitCommand = [
 		'sbatch',
 		'-J{0}'.format('Antigen'),
 		'-D{0}'.format(resultsDir),
 		'--array={0}'.format(','.join(str(i) for i in JOBS)),
 		sbatchFilename
 	]
	
	# Print command to terminal
	process = subprocess.Popen(['echo'] + submitCommand)
	process.wait()
	
	# Construct environment variables
	env = dict(os.environ)
	env['ANTIGEN_ROOT'] = rootDir
	
	# Actually run command
	if not DRY:
		process = subprocess.Popen(submitCommand, env=env)
		process.wait()
        
if __name__ == '__main__':
    sweepDir = os.path.abspath(os.path.dirname(__file__))
    rootDir = os.path.abspath(os.path.join(sweepDir, '..'))
    resultsDir = os.path.join(sweepDir, 'results')
    submitJob(rootDir, sweepDir, resultsDir)
