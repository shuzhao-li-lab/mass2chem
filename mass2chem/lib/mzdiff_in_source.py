'''

will use https://github.com/stanstrup/commonMZ for now, which come mostly from Keller et al, 2008.

Could consider the tables from CAMERA later.

Two sets of tables for pos and neg ESI.

'''
import os
from collections import namedtuple
from itertools import combinations

import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity

from pyopenms import MSExperiment, MzMLFile

# Define major data strctures

Scan = namedtuple('Scan', [
    'scan_number', 'retention_time', 'mz_values', 'intensity_values'
    ])

class Sample: # modified from asari
    '''
    A sample or injection, corresponding to one raw data file, either from positive or negative ionization.
    '''
    def __init__(self, experiment, mode, input_file=''):
        self.input_file = input_file
        self.experiment = experiment                        # parent Experiment instance
        self.name = os.path.basename(input_file)            # must be unique
        self.mode = mode            # ESI mode, could switch within expt
        
        # list of scans, each as one spectrum collected at a specific elution time
        self.scans = []
        
        # these two lists are in mached order; will make immutable tuples after reading data
        self.retention_index = []     # this corresponds to scan number
        self.retention_time = []
        


def get_pairwise_diffs(L):
    # limit to diff btw (1,100) for efficiency. E.g. scan #1 234955 -> 
    d = [abs(x-y) for x,y in combinations(L,2)]
    return [x for x in d if 1<x<100]

def read_sample(input_file, expt_='ecoli', mode='pos'): 
    SM = Sample(expt_, mode, input_file) 
    exp = MSExperiment()
    MzMLFile().load(input_file, exp)
    print(input_file, exp.getNrSpectra())
    scans = []
    ii = 0
    for sp in exp:
        ms_level = sp.getMSLevel()
        if ms_level == 1:    # ONLY dealig with MS1 here
            rt = sp.getRT()
            mz, intensity = sp.get_peaks()
            # use int and tuple to save space when storage is considered
            intensity = [int(x) for x in intensity]
            # Scan defined as 'scan_number', 'retention_time', 'mz_values', 'intensity_values'
            scans.append(Scan(ii, rt, tuple(mz), tuple(intensity)))
            # scans.append(Scan(ii, rt, mz, intensity))
        ii += 1
        
    SM.scans = scans
    SM.retention_index = tuple([sc.scan_number for sc in scans])
    SM.retention_time = tuple([sc.retention_time for sc in scans])
    return SM








