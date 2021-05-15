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

maxJobs = 200

NUM_RUNS = 20
VACCINATIONRATES = np.array([0,0.01,0.05, 0.1, 0.2, 0.3])
totalJobs = NUM_RUNS*len(VACCINATIONRATES)

N_PER_JOB = totalJobs/maxJobs+1
if totalJobs <= maxJobs:
    N_PER_JOB = 1
    maxJobs = totalJobs

totalN = 50000000
birthRate = 0.000091
Reff = 1.0
R0 = 1.8
nu = 0.2
beta = R0 * (nu+birthRate)
meanStep = 0.6
sdStep = 0.3

vaccinate = True

muPhenotype = 0
useReferenceStrains = True

vaccineImmuneBreadth = 1
vaccineLag = 300
deployDay = vaccineLag
vaccineWindow = 120
vaccinateConstantFraction = True
startAtEquilibriumInfected = True
startAtEquilibriumImmune = True
fractionNeverVaccinated = 0.33
fractionRepeatVaccinations = 0.8

referenceIDs = [1,2,5,6,9,11,15,16,17,26,29,30,31,32,33,36,37,38,39,41,43,44,48,49,50,51,53,54,55,58,59,62,63,64,65,68,69,72,73,74,75,80,82,84,85,86,87,90,91,92,93,94,96,97,98,100,101,102,106,107,110,118,119,120,122,124,125,126,127,128,130,131,134,135,136,137,138,141,143,146,147,148,149,150,151,153,154,155,157,158,159,161,162,163,164,166,167,168,172,173,178,183,184,185,187,190,191,195,197,199,202,203,204,205,206,208,209,211,217,218,219,225,227,232,233,234,235,236,240,242,243,245,246,247,248,250,253,256,258,259,263,265,267,270,275,276,277,278,279,280,281,282,283,287,289,290,291,292,295,296,297,298,299,300,304,306,307,310,311,312,313,319,321,325,326,328,331,332,333,334,335,337,338,341,342,343,344,346,348,352,353,354,359,360,361,362,363,364,365,367,368,370,371,372,376,380,383,388,390,393,395,396,397,399,402,404,405,406,409,410,413,414,415,417,418,421,422,423,424,425,426,428,430,431,432,433,435,436,438,440,441,442,443,444,445,447,448,449,452,455,456,461,464,465,469,472,473,474,478,479,480,482,490,491,497]
def generateJobs():    
    jobNum = 0
    
    for vaccinationRate in VACCINATIONRATES:
        for runNum in range(NUM_RUNS):
            # This is the name that SLURM uses to identify the job
            jobName = jobNum
            # This is the subdirectory inside 'results' used to run the job
            jobSubdir = '{}'.format(jobName)
            
            referenceFile = '../../../../make_reference/vaccine/results/' + str(np.random.choice(referenceIDs,1)[0]) + '/reference.sqlite'

            paramDict = OrderedDict([
                ('referenceFile', referenceFile),
                ('useReferenceStrains', useReferenceStrains),
                ('muPhenotype', [muPhenotype]),
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
    sbatchFilename = os.path.join(sweepDir, 'vaccine_longitudinal.sbatch')
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
