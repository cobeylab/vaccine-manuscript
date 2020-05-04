#!/usr/bin/env python

# jobs are indexed staring from 0

import csv
import os
import json
import numpy as np
import subprocess
from collections import OrderedDict

###############################################################################

### Modify these values and the generateJobs() function ###

# Set DRY = False to actually submit jobs
DRY = False

maxJobs = 500
NUM_RATES = 40
R0S = np.linspace(1.2,3,10)

NUM_RUNS = 40
totalJobs = NUM_RUNS*len(R0S)*NUM_RATES
if totalJobs < maxJobs:
    maxJobs = totalJobs
N_PER_JOB = totalJobs/maxJobs+1

vaccinate = [True]

vaccineImmuneBreadth = 1
vaccineLag = 300
deployDay = vaccineLag
vaccineWindow = 120

gamma = 0.2
birthRate = .000091

def generateJobs():    
    jobNum = 0
    
    for R0 in R0S:  #unvaccinated R0 to be consistent with other tests
        beta = (gamma+birthRate)*R0
        thresholdRate = -(gamma+2*birthRate)/2 + (gamma**2 + 4*birthRate*beta)**.5/2
        rateLB = thresholdRate - 0.006/365.0
        rateUB = thresholdRate + 0.004/365.0
        VACCINATIONRATES = np.linspace(rateLB, rateUB, NUM_RATES)
        for vaccinationRate in VACCINATIONRATES:
            #calibrate equilibria to beta instead of vaccine-free R0
            vaccinatedR0 = beta/(gamma+birthRate + vaccinationRate)
            eqS = 1/vaccinatedR0
            eqI = birthRate/beta*(vaccinatedR0-1)  - vaccinationRate/beta
            eqR = 1- eqS - eqI

            if eqI <= 0:
                eqR = vaccinationRate/(birthRate + vaccinationRate)
                eqI = 1/5e7
                eqS = 1-eqR-eqI
                
            initialPrR = eqR
            
            for runNum in range(NUM_RUNS):
                # This is the name that SLURM uses to identify the job
                jobName = jobNum

                # This is the subdirectory inside 'results' used to run the job
                jobSubdir = '{}'.format(jobName)

                paramDict = OrderedDict([
                    ('vaccinationRate', [vaccinationRate]),
                    ('initialPrR', initialPrR),
                    ('initialIs', [int(eqI*50000000)]),
                    ('initialRs', [int(initialPrR*50000000)]),
                    ('vaccinationRate', [vaccinationRate]),
                    ('startAtEquilibriumInfected', [False]),
                    ('startAtEquilibriumImmune', [False]),
                    ('vaccineImmuneBreadth',vaccineImmuneBreadth),
                    ('vaccineLag', vaccineLag),
                    ('vaccineWindow', vaccineWindow),
                    ('deployDay', deployDay),
                    ('R0',R0),
                    ('beta',beta)
                ])
                yield (jobName, jobSubdir, paramDict)
                jobNum += 1


###############################################################################

def writeParameters(resultsDir, jobSubdir, paramDict):
    jobDir = os.path.join(resultsDir, jobSubdir)
    os.makedirs(jobDir)
    paramsFilename = os.path.join(jobDir, 'parameters.json')
    with open(paramsFilename, 'w') as paramsFile:
        json.dump(paramDict, paramsFile, indent=2)
        paramsFile.write('\n')

def submitJob(rootDir, sweepDir, resultsDir):
    # Construct SLURM command
    sbatchFilename = os.path.join(sweepDir, 'vaccine.sbatch')
    submitCommand = [
        'sbatch',
        '-J{0}'.format('Antigen'),
        '-D{0}'.format(resultsDir),
        #'--array=0-1',
        '--array=0-{0}'.format(maxJobs-1),
        sbatchFilename
    ]
    
    # Print command to terminal
    process = subprocess.Popen(['echo'] + submitCommand)
    process.wait()
    
    # Construct environment variables
    env = dict(os.environ)
    env['ANTIGEN_ROOT'] = rootDir
    env['N_PER_JOB'] = str(N_PER_JOB)
    
    # Actually run command
    if not DRY:
        process = subprocess.Popen(submitCommand, env=env)
        process.wait()

def loadUncommentedJsonString(filename):
    lines = list()
    with open(filename) as jsonFile:
        for line in jsonFile:
            lines.append(line.rstrip('\n').split('//')[0])
    return '\n'.join(lines)

if __name__ == '__main__':
    sweepDir = os.path.abspath(os.path.join(os.path.dirname(__file__)))
    rootDir = os.path.abspath(os.path.join(sweepDir, '..'))
    resultsDir = os.path.join(sweepDir, 'results')
    
    cParamsFilename = os.path.join(sweepDir, 'vaccine_parameters.json')
    cParamsJson = loadUncommentedJsonString(cParamsFilename)
    cParamsDict = json.loads(cParamsJson, object_pairs_hook=OrderedDict)
    
    for jobName, jobSubdir, paramDict in generateJobs():
        allParamsDict = OrderedDict()
        allParamsDict.update(cParamsDict)
        allParamsDict.update(paramDict)
        writeParameters(resultsDir, jobSubdir, allParamsDict)
        
    submitJob(rootDir, sweepDir, resultsDir)
