#!/usr/bin/env python

# jobs are indexed staring from 0
import sys
import csv
import os
import json
import numpy as np
import subprocess
from collections import OrderedDict
import pandas as pd
import sqlite3

PARAMS = ['randomSeed', 'muPhenotype', 'meanStep', 'sdStep', 'initialTraitA','beta','smithConversion','endDay','vaccineWindow','vaccineSD','vaccinationRate','vaccineImmuneBreadth','vaccineLag']
PARAM_TYPES = {
    'randomSeed' : 'REAL',
    'muPhenotype' : 'REAL',
    'meanStep' : 'REAL',
    'sdStep' : 'REAL',
    'initialTraitA' : 'REAL',
    'beta' : 'REAL',
    'smithConversion' : 'REAL',
    'endDay' : 'REAL',
    'vaccineWindow' : 'REAL',
    'vaccineSD' : 'REAL',
    'vaccinationRate' : 'REAL',
    'vaccineImmuneBreadth' : 'REAL',
    'vaccineLag' : 'REAL',
}

VACPARAMS = ['beta','nu','vaccineSD','vaccinationRate','vaccineImmuneBreadth','vaccineLag','deployDay','vaccineWindow']
VACPARAM_TYPES = {
    'beta' : 'REAL',
    'nu' : 'REAL',
    'vaccineSD' : 'REAL',
    'vaccinationRate' : 'REAL',
    'vaccineImmuneBreadth' : 'REAL',
    'vaccineLag' : 'REAL',
    'deployDay' : 'REAL',
    'vaccineWindow' : 'REAL',
}

def loadUncommentedJsonString(filename):
    lines = list()
    with open(filename) as jsonFile:
        for line in jsonFile:
            lines.append(line.rstrip('\n').split('//')[0])
    return '\n'.join(lines)

def insert(db, table, values):
    return db.execute(
        '''INSERT INTO {0} VALUES ({1})'''.format(
            table,
            ','.join(['?'] * len(values))
        ),
        values
    )
    
def loadJob(comboDb, resultsDir, runId, jobDir):
    # Load dictionary of parameter values from parameters.json file
    paramsJson = loadUncommentedJsonString(os.path.join(jobDir, 'parameters_out.json'))
    paramsDict = json.loads(paramsJson, object_pairs_hook=OrderedDict)
    paramVals = paramsDict

    insert(comboDb, 'parameters', [runId] + [paramVals[paramName] if type(paramVals[paramName]) != list else paramVals[paramName][0] for paramName in PARAMS])
    comboDb.commit()
    
    outDb = sqlite3.connect(os.path.join(jobDir,'output.sqlite'))
    
    timeseries = pd.read_sql('SELECT date, totalN, totalCases FROM timeseries', con=outDb)
    meanAnnualIncidence = np.sum(timeseries['totalCases']) / max(timeseries['date']) / np.mean(timeseries['totalN'])
    cumulativeIncidence = np.sum(timeseries['totalCases'])
    lastDate = max(timeseries['date'])
    
    
    extinct = outDb.execute('SELECT went_extinct FROM status').next()[0]
    tmrcaLimit = outDb.execute('SELECT reached_tmrca_limit FROM status').next()[0]
    fluLike = 1-max(extinct,tmrcaLimit)
    insert(comboDb, 'status',[runId, extinct,tmrcaLimit,fluLike])
        
    outDb.close()
    comboDb.commit()

def loadJobs(comboDb, sweepDir, resultsDir):

    subDirs = os.listdir(resultsDir)
    runId = 0
    for subDir in subDirs:
        jobDir = os.path.join(resultsDir, subDir)
        runId = subDir
        #print(str(runId) + ' : ' + str(subDirs.index(subDir)))
        if not os.path.isfile(os.path.join(jobDir, 'parameters_out.json')):
            print 'missing parameters_out.json ' + str(runId)
        if not os.path.isfile(os.path.join(jobDir, 'output.sqlite')):
            print 'missing output.sqlite ' + str(runId)
        if os.path.isdir(jobDir) & os.path.isfile(os.path.join(jobDir, 'parameters_out.json')) & os.path.isfile(os.path.join(jobDir, 'output.sqlite')):
            try:
                loadJob(comboDb, resultsDir, runId, jobDir)
            except:
                e = sys.exc_info()[0]
                print "ERROR " + str(runId)
                #print ("<p>Error: %s</p>" % e)
                

if __name__ == "__main__":
    sweepDir = os.path.abspath(os.path.dirname(__file__))
    resultsDir = os.path.join(sweepDir, 'results')
    
    
    comboDb = sqlite3.connect(os.path.join(sweepDir, 'results.sqlite'))
    
    # Index of runs with parameters
    comboDb.execute(
        "CREATE TABLE parameters (runId INTEGER, {0})".format(
            ', '.join([paramName + ' ' + PARAM_TYPES[paramName] for paramName in PARAMS])
        )
    )
    comboDb.execute("CREATE TABLE cumulativeIncidence (runId INTEGER, cumulativeIncidence REAL);")
    comboDb.execute("CREATE TABLE cumulativeDrift (runId INTEGER, culumativeDrift REAL);")
    comboDb.execute("CREATE TABLE status (runId INTEGER, extinct INTEGER, excessDiversity INTEGER, fluLike INTEGER);")
    comboDb.execute("CREATE TABLE tipFluxRate (runId INTEGER, tipFluxRate REAL);")
    comboDb.execute("CREATE TABLE meanFluxRate (runId INTEGER, meanFluxRate REAL);")
    comboDb.execute("CREATE TABLE meanAnnualIncidence (runId INTEGER, meanAnnualIncidence REAL);")
    comboDb.execute("CREATE TABLE pooled_results (runId INTEGER, cumulativeIncidence REAL, cumulativeDrift REAL, meanAnnualIncidence REAL, meanFluxRate REAL, tipFluxRate REAL, lastDate REAL, meanTMRCA REAL, {0}, extinct REAL, tmrcaLimit REAL, fluLike REAL)".format(
            ', '.join([paramName + ' ' + VACPARAM_TYPES[paramName] for paramName in VACPARAMS])
        )
    )
    
    comboDb.commit()    
    loadJobs(comboDb, sweepDir, resultsDir)

    comboDb.close()
        
