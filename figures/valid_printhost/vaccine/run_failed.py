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
        '--array=99',
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
