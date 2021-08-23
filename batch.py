#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 14:29:22 2021

@author: ja17375

batch.py enables batch processing of shear-wave splitting using a Python wrapper
"""

from multiprocessing import Pool, current_process
from functools import partial
import contextlib
import pandas as pd 
import obspy
from wrapper import Wrapper

SHEBA_EXEC = '/Users/ja17375/Ext_Programs/bin'

def measure_event(file, rundir, phase):
    '''
    Measures shear-wave splitting for a given event
    '''
    st = obspy.read(file)
    Measure = Wrapper(st, phase, rundir=rundir)
    Measure.preprocess()
    outfile = file.split('.')[0].split('/')[-1] # selects just filename with not extentsion or preceding path
    result = Measure.measure_splitting(outfile, SHEBA_EXEC)

def serial_process(filelist, rundir, phase):
    '''
    Measures shear-wave splitting for all files in filelist. Iterates over filelist.
    '''

    for file in filelist:
        measure_event(file, rundir, phase)

def pool_process(filelist, rundir, phase, ncores=2):
    '''
    Measures shear-wave splitting using parallel processing (Pool)
    '''
    runner = partial(measure_event,rundir,phase)
    print("Parallel job. {} cores requested".format(ncores))
    with contextlib.closing( Pool(processes = ncores) ) as pool:
    #           Iterate over stations in the station list.
        pool.map(runner,filelist)