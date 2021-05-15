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
VACCINATIONRATES = np.linspace(0,0.5,11)[1:11]
VACCINEIMMUNEBREADTHS = np.array([0.05, 0.1, 0.2, 0.3, 0.5 , 0.7])
VACCINEIMMUNEBREADTHS = np.array([1])

NUM_RUNS = 200
totalJobs = NUM_RUNS*len(VACCINATIONRATES)*len(VACCINEIMMUNEBREADTHS)
if totalJobs < maxJobs:
    maxJobs = totalJobs
N_PER_JOB = totalJobs/maxJobs+1

totalN = 50000000
birthRate = 0.000091
Reff = 1.0
R0 = 1.8
nu = 0.2
beta = R0 * (nu+birthRate)
meanStep = 0.6
sdStep = 0.3

vaccinate = True

#vaccineImmuneBreadth = 1
vaccineLag = 0
deployDay = vaccineLag
vaccineWindow = 120
vaccinateConstantFraction = True
startAtEquilibriumInfected = True
startAtEquilibriumImmune = True
fractionNeverVaccinated = 0.33
fractionRepeatVaccinations = 0.8



def generateJobs():
    jobNum = 0
    for vaccineImmuneBreadth in VACCINEIMMUNEBREADTHS:
        for vaccinationRate in VACCINATIONRATES:
    #         initialIs = (birthRate/beta * (R0 * (1-vaccinationRate) - 1 )) * totalN
    #         initialIs = int(initialIs)
    #         if initialIs < 0:
    #             initialIs = 10

            for runNum in range(NUM_RUNS):
                # This is the name that SLURM uses to identify the job
                jobName = jobNum
                # This is the subdirectory inside 'results' used to run the job
                jobSubdir = '{}'.format(jobName)

                paramDict = OrderedDict([
                    ('vaccinate', vaccinate),
                    ('vaccinationRate', [vaccinationRate]),
                    ('vaccineImmuneBreadth',vaccineImmuneBreadth),
                    ('vaccineLag', vaccineLag),
                    ('vaccineWindow', vaccineWindow),
                    ('deployDay', deployDay),
                    ('meanStep',meanStep),
                    ('sdStep',sdStep),
                    ('R0',R0),
                    ('beta',beta),
                    ('vaccinateConstantFraction', vaccinateConstantFraction),
                    ('startAtEquilibriumInfected', [startAtEquilibriumInfected]),
                    ('startAtEquilibriumImmune', [startAtEquilibriumImmune]),
                    ('fractionNeverVaccinated', fractionNeverVaccinated),
                    ('fractionRepeatVaccinations', fractionRepeatVaccinations)
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
