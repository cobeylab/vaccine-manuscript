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

def write_timeseries(job_id, jobDir, out_db):
    fileName = os.path.join(jobDir,'out.timeseries')
    if(os.path.isfile(fileName)):
		with open(os.path.join(jobDir,'out.timeseries'), 'rU') as ts_file:
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

def write_branches(job_id, jobDir, out_db):
    fileName = os.path.join(jobDir,'out.branches')
    if(os.path.exists(fileName)):
		with open(os.path.join(jobDir,'out.branches'), 'rU') as ts_file:
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

def write_tips(job_id, jobDir, out_db):
	fileName = os.path.join(jobDir,'out.tips')
	if(os.path.exists(fileName)):
		with open(os.path.join(jobDir,'out.tips'), 'rU') as file:
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

def write_trees(job_id, jobDir, out_db):
	fileName = os.path.join(jobDir,'out.trees')
	if(os.path.exists(fileName)):
		with open(os.path.join(jobDir,'out.trees'), 'rU') as file:
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
    went_extinct = os.path.exists(os.path.join(jobDir,'out.extinct'))
    reached_tmrca_limit = os.path.exists(os.path.join(jobDir,'out.tmrcaLimit'))
    insert(out_db, 'status', [job_id, went_extinct, reached_tmrca_limit])

def loadJobs(sweepDir, resultsDir, out_db):
	# Extract parameters names
	subDirs = os.listdir(resultsDir)
	runId = 0
	for subDir in subDirs:
		jobDir = os.path.join(resultsDir, subDir)
		runId = subDir
		if os.path.isdir(jobDir):
			write_status(subDir,jobDir,out_db)
			write_timeseries(subDir,jobDir,out_db)
			write_branches(subDir,jobDir,out_db)
			write_tips(subDir,jobDir,out_db)
			#write_trees(subDir,jobDir,out_db) #trees are too large for field size
			out_db.commit()
			

if __name__ == '__main__':
	sweepDir = os.path.abspath(os.path.dirname(__file__))
	resultsDir = os.path.join(sweepDir, 'results')
	
	out_db = sqlite3.connect(os.path.join(sweepDir, 'output.sqlite'))
	
	loadJobs(sweepDir, resultsDir, out_db)

	out_db.close()
