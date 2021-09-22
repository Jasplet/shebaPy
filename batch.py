#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 14:29:22 2021

@author: ja17375

batch.py enables batch processing of shear-wave splitting using a Python wrapper
"""

#from multiprocessing import Pool, current_process
# from functools import partial
# import contextlib
import pandas as pd
import obspy
from wrapper import Wrapper

SHEBA_EXEC = '/Users/ja17375/Ext_Programs/bin'

def measure_event(file, rundir, phase, window=False):
    '''
    Measures shear-wave splitting for a given event
    '''
    input_data = obspy.read(file)
    sheba_wrapper = Wrapper(input_data, phase, rundir=rundir)
    sheba_wrapper.preprocess()
    outfile = file.split('.')[0].split('/')[-1]
    # selects just filename with not extentsion or preceding path
    splitting_result = sheba_wrapper.measure_splitting(outfile, SHEBA_EXEC, window)
    return splitting_result

def serial_process(filelist, rundir, phases, window=False):
    '''
    Measures shear-wave splitting for all files in filelist. Iterates over filelist.
    '''
    results = []
        
    for i, file in enumerate(filelist):
            
        result = measure_event(file, rundir, phases[i], window)
        results.append(result)

    result_df = pd.DataFrame(results)
    return result_df

# def pool_process(filelist, rundir, phase, ncores=2):
#     '''
#     Measures shear-wave splitting using parallel processing (Pool)
#     '''
#     runner = partial(measure_event,rundir,phase)
#     print("Parallel job. {} cores requested".format(ncores))
#     with contextlib.closing( Pool(processes = ncores) ) as pool:
#         pool.map(runner,filelist)      