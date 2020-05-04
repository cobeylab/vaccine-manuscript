#!/usr/bin/env python

import os
import json
import sqlite3
import csv
import numpy as np
import subprocess
from collections import OrderedDict
from datetime import datetime
import random
import time
import argparse
import sys

MATHOUTS = ['antigenicFlux', 
			'antigenicLagAG1', 
			'meanAnnualDemeIncidence', 
			'meanAnnualIncidence', 
			'meanDemeDiversity', 
			'meanDemeNetau',
			'meanDemeTMRCA',
			'meanDiversity', 
			'meanFluxRate',
			'meanNeTau',
			'meanTMRCA',
			'mutationSizeFlux',
			'prevalenceSpatial',
			'sideBranchMigRate',
			'sideMigrationRateMatrix',
			'spatialMutationRate',
			'spatialMutationSize',
			'trunkMigRate',
			'trunkMigrationRateMatrix',
			'trunkProportions']

def isfloat(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

def insert(db, table, values):
    return db.execute(
        '''INSERT INTO {0} VALUES ({1})'''.format(
            table,
            ','.join(['?'] * len(values))
        ),
        values
    )

def write_summary(job_id, jobDir, out_db):
	fileName = os.path.join(jobDir, 'out.summary')
	if(os.path.isfile(fileName)):
		with open(fileName, 'rU') as summaryFile:
			cw = csv.reader(summaryFile, delimiter='\t')
			cw.next()
			out_db.execute(
				'''CREATE TABLE IF NOT EXISTS summary (
					runId INTEGER, key TEXT, value TEXT
				)'''
			)
			for key, value in cw:
				insert(out_db,'summary',[job_id, key, value])
		os.remove(fileName)


def write_timeseries(job_id, jobDir, out_db):
    fileName = os.path.join(jobDir,'out.timeseries')
    if(os.path.isfile(fileName)):
		with open(fileName, 'rU') as ts_file:
			cr = csv.reader(ts_file, delimiter='\t')
			header = cr.next()
			out_db.execute(
				'''CREATE TABLE IF NOT EXISTS timeseries (
					job_id INTEGER,
					{0}
				)'''.format(', '.join(['{0} REAL'.format(x) for x in header]))
			)
			for row in cr:
				insert(out_db, 'timeseries', [job_id] + [float(x) for x in row])
		os.remove(fileName)

def write_branches(job_id, jobDir, out_db):
    fileName = os.path.join(jobDir,'out.branches')
    if(os.path.exists(fileName)):
		with open(fileName, 'rU') as ts_file:
			cr = csv.reader(ts_file, delimiter='\t')
			header = ['name','year','trunk','tip','mark','location','layout','ag1','ag2',
				'parentname','parentyear','parenttrunk','parenttip','parentmark','parentlocation','parentlayout','parentag1','parentag2',
				'coverage']
			out_db.execute(
				'''CREATE TABLE IF NOT EXISTS branches (
					job_id INTEGER,
					{0}
				)'''.format(', '.join(['{0} REAL'.format(x) for x in header]))
			)
			for row in cr:
				row = row[0].split(',') + row[1].split(',') + row[2].split(',')
				row = [x.replace('{','').replace('}','') for x in row]
				insert(out_db, 'branches', [job_id] + [float(x) if isfloat(x) else x for x in row])
		os.remove(fileName)

def write_tips(job_id, jobDir, out_db):
	fileName = os.path.join(jobDir,'out.tips')
	if(os.path.exists(fileName)):
		with open(fileName, 'rU') as file:
			cr = csv.reader(file, delimiter=',')
			header = cr.next()
			out_db.execute(
				'''CREATE TABLE IF NOT EXISTS tips (
					job_id INTEGER,
					{0}
				)'''.format(', '.join(['{0} REAL'.format(x) for x in header]))
			)
			for row in cr:
				insert(out_db, 'tips', [job_id] + [float(x) if isfloat(x) else x for x in row])
		os.remove(fileName)


def write_trees(job_id, jobDir, out_db):
	fileName = os.path.join(jobDir,'out.trees')
	if(os.path.exists(fileName)):
		with open(fileName, 'rU') as file:
			cr = csv.reader(file, delimiter=' ')
			header = ['tree']
			out_db.execute(
				'''CREATE TABLE IF NOT EXISTS trees (
					job_id INTEGER,
					{0}
				)'''.format(', '.join(['{0} TEXT'.format(x) for x in header]))
			)
			for row in cr:
				insert(out_db, 'trees', [job_id] + [x for x in row])

def write_status(job_id, jobDir, out_db):
    out_db.execute('CREATE TABLE IF NOT EXISTS status (job_id INTEGER, went_extinct INTEGER, reached_tmrca_limit INTEGER)')
    extinctFile = os.path.join(jobDir,'out.extinct')
    tmrcaFile = os.path.join(jobDir,'out.tmrcaLimit')
    went_extinct = os.path.exists(extinctFile)
    reached_tmrca_limit = os.path.exists(tmrcaFile)
    insert(out_db, 'status', [job_id, went_extinct, reached_tmrca_limit])
    if(went_extinct):
    	os.remove(extinctFile)
    if(reached_tmrca_limit):
    	os.remove(tmrcaFile)

def loadCsv(mathOut,comboDb,jobDir):
	
	fileName = os.path.join(jobDir, mathOut+'.csv')
	if(os.path.exists(fileName)):
		with open(fileName) as mathFile:
			cw = csv.reader(mathFile, delimiter=',')
			ncols = len(next(cw))
			mathFile.seek(0)
			for row in cw:
				length = len(row)
				if length > ncols:
					ncols = length
			mathFile.seek(0)
			HEADERS = [str(n) for n in list(xrange(ncols))]
			TYPES = ['TEXT']*ncols
			cols = ','.join(['"%s" %s' % (header, type) for (header,type) in zip(HEADERS, TYPES)])
			try:
				comboDb.execute(
					"CREATE TABLE {0} ({1});".format(mathOut, cols)
				)
				for line in cw:
					#print "INSERT INTO {0} VALUES ({1})".format(mathOut,', '.join([str(n) for n in line]))
					try:
						comboDb.execute(
							"INSERT INTO {0} VALUES ({1})".format(mathOut,','.join('?'*ncols)),[str(n) for n in line]
						)
					except:
						pass
			except:
				pass
		comboDb.commit()
		os.remove(fileName)

if __name__ == '__main__':
	sweepDir = os.path.abspath(os.path.dirname(__file__))
	resultsDir = os.path.join(sweepDir, 'results')
	subDir = sys.argv[1]
	jobDir = os.path.join(resultsDir, subDir)
	if(os.path.isfile(os.path.join(jobDir,'out.bashscreen'))):
		os.remove(os.path.join(jobDir,'out.bashscreen'))
	if(os.path.isfile(os.path.join(jobDir,'out.range'))):
		os.remove(os.path.join(jobDir,'out.range'))
	
	#if os.path.isfile(os.path.join(jobDir, 'out.branches')):
	#	comboDb = sqlite3.connect(os.path.join(jobDir, 'math.sqlite'))
	#	for mathOut in MATHOUTS:
	#		#cleanCsv(mathOut,jobDir)
	#		loadCsv(mathOut,comboDb,jobDir)
	#	comboDb.commit()
	#	comboDb.close()

	out_db = sqlite3.connect(os.path.join(jobDir, 'output.sqlite'))
	write_summary(subDir,jobDir,out_db)
	write_status(subDir,jobDir,out_db)
	write_timeseries(subDir,jobDir,out_db)
	write_branches(subDir,jobDir,out_db)
	write_tips(subDir,jobDir,out_db)
	out_db.commit()
	out_db.close()

