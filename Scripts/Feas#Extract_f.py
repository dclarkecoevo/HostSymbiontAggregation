


from __future__ import division
import agg_fxns as agg
import pandas as pd
import macroeco.compare as comp
import scipy.stats as stats
import pypartitions as pyp
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Define file locations.
outputPath = 'OutputFeas_Set_FemaleNew.csv'
inputPath = 'FINAL_TT1_NEW12.csv'

if os.path.exists(outputPath):
    os.remove(outputPath)

# Initialize data.
data = pd.read_csv(inputPath)
print('Loaded ' + inputPath + '.')

# This dict contains data subsetted to 'nsiteyear100' with each column existing as a key 
# with an array as the value in the key/value pair.
# Example: { "caroni:arima:lower:1:2003" : { "gyro": [1,1,1,1], "sex": ["f","f","j","m"], "Length" : [12, 11, 10, 9]} }
holding = {}

print('Building holding dictionary...')

# Iterate through all rows 
for i, row in data.iterrows():
    if row.sex == 'f':

        # Check if the site is being held, if not create a key for it in the dict.
        if not holding.has_key(row.nsiteyear):
            holding[row.nsiteyear] = {}

        # For each column of that site, push the value of the column to the corresponding dict.
        for column in row.index:
            if not holding[row.nsiteyear].has_key(column):
                 holding[row.nsiteyear][column] = []
            holding[row.nsiteyear][column].append(row[column])

print('Generating data frames...\n')

# Define debug vars.
count = 0
length = len(holding)

# For each site, run the calc and print the results.
for key in holding:

    # Debug - Can be deleted.
    count = count+1
    percent = float(count) * 100 / length
    bar = '|' * int(percent/100 * 20 - 1) + '|'
    spaces = ' ' * (20 - len(bar))
    sys.stdout.write('Generation progress - [%s%s] %d%%\r' % (bar, spaces, percent))
    if count != length:
        sys.stdout.flush()
    else:
         sys.stdout.write('\n')
    
    # Run calculcations.
    emp = ','.join(map(str, np.sort(holding[key]['gyro'])))
    H = len(holding[key]['gyro'])
    P = np.sum(holding[key]['gyro'])
    feas = ','.join(map(str, agg.feasible_mixture([(P, H)],samples=1000, center="median")[1]))
    ranks = ','.join(map(str, np.arange(1, H + 1)))

    # Generate output.
    output = pd.DataFrame({'feas' : feas,
                           'emp' : emp,
                           'ranks' : ranks,
                           'column' : key},
             index = [0],
             columns=['feas','emp','ranks','column'])
    output.to_csv(outputPath, index=False, mode='a', line_terminator='\n', header=not os.path.exists(outputPath)) 

print('Generation complete.')
