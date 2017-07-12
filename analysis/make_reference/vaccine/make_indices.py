#!/usr/bin/env python

import os
import sys
import csv
import sqlite3
import json
import numpy as np
import subprocess
from collections import OrderedDict

if __name__ == "__main__":
    sweepDir = os.path.abspath(os.path.dirname(__file__))
    resultsDir = os.path.join(sweepDir, 'results')
    jobDirs = range(175,200)
    jobDirs = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,50,51,52,54,55,56,57,58,59,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,81,84,85,86,87,88,89,90,91,92,93,94,95,96,97,99,100,102,105,106,108,112,113,117,118,119,120,121,129,132,137,139,175,401,402,403,404,405,406,407,409,410,411,412,414,415,416,419,422,424]
    jobDirs = [173,174]
    jobDirs = range(0,500)
 
    for jobId in jobDirs:
        print(jobId)
        jobDir = os.path.join(resultsDir, str(jobId))
        comboDb = sqlite3.connect(os.path.join(jobDir, 'reference.sqlite'))
        
        try:
            # Index of runs with parameters
            comboDb.execute(
                "CREATE INDEX `day` ON `phenotypes` (`day` )"
            )
        except:
            pass
    
        comboDb.commit()
    
        comboDb.close()
    
