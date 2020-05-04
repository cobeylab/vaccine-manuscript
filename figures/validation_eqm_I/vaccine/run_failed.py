#!/usr/bin/env python

import csv
import os
import json
import numpy as np
import subprocess
from collections import OrderedDict

DRY = False

def submitJob(rootDir, sweepDir, resultsDir):
    # Construct SLURM command
    sbatchFilename = os.path.join(sweepDir, 'run_failed.sbatch')
    
    submitCommand = [
        'sbatch',
        '-J{0}'.format('Antigen'),
        '-D{0}'.format(resultsDir),
        '--array=80,88,199,107,83,187,171,176,92,114,110,190,87,169,161,195,191,4,1,178,112,0,104,100,180,86,179,81,89,84,101,185,95,90,186,98,162',
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
