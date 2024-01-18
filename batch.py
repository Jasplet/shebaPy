#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 14:29:22 2021

@author: ja17375

batch.py enables batch processing of shear-wave splitting using a Python wrapper
"""

import pandas as pd
import obspy
import os
#Imports for multi-processing
#from multiprocessing import Pool, current_process
# from functools import partial
# import contextlib
#Local import within package
from wrapper import Wrapper
from wrapper import collate_result
SHEBA_EXEC = '/Users/eart0593/Ext_Programs/bin'

def measure_event(file, rundir, phase, window=False, nwind=10, tlag_max=4.0, debug=False, c1=0.01, c2=0.5, trim=True):
    '''
    Measures shear-wave splitting for a given event
    
    Parameters
    ----------
    file : str 
        name of SAC files containing waveform data to measure
    rundir : str
        directory in which Sheba will run to make the measurement
    phase : str
        shear-wave phase to meausre
    window : boolean, optional, default=False
        switch to enabale manual interactive windowing (if True)
    
    Returns
    ----------
    splitting_result : dict
        shear-wave splitting result and metadata
    '''
    input_data = obspy.read(file)
    sheba_wrapper = Wrapper(input_data, phase, rundir=rundir)
    sheba_wrapper.preprocess(c1, c2, trim)
    fstem = file.split('/')[-1]
    outfile = '.'.join(fstem.split('.')[0:-1])

    # selects just filename with not extentsion or preceding path
    splitting_result = sheba_wrapper.measure_splitting(outfile, SHEBA_EXEC, window, nwind, tlag_max, debug)
    return splitting_result

def serial_process(filelist, rundir, phases, window=False, nwind=10, tlag_max=4.0, debug=False, c1=0.01, c2=0.5, trim=True, remeasure=False):
    '''
    Measures shear-wave splitting for all files in filelist. Iterates over filelist.
    
    Parameters
    ----------
    filelist : list 
        list of SAC files (as strings) containing waveform data to measure
    rundir : str
        directory in which Sheba will run to make the measurement
    phases : list
        list of shear-wave phases to meausre for each file
    window : boolean, optional, default=False
        switch to enabale manual interactive windowing (if True)
    
    Returns
    ----------
    result_df : DataFrame
        shear-wave splitting results and metadata collated into a pandas Dataframe  
    '''
    results = []


    for i, file in enumerate(filelist):
        res_stem = file.split('/')[-1]
        res_id = res_stem.split('.')[0]
        if os.path.isfile(f'{rundir}/{res_id}_sheba_result.nc'):
            if remeasure:
                try:
                    result = measure_event(file, rundir, phases[i], window, nwind, tlag_max, debug, c1, c2, trim)
                    if result is not None:
                        #Only append full Dicts!
                        results.append(result)
                except ValueError:
                    print(ValueError)
            else: 
                print('Skipping event - already measured')
                result = collate_result(rundir, res_id)
                results.append(result)
        else:
            try:
                print(file)
                result = measure_event(file, rundir, phases[i], window, nwind, tlag_max, debug, c1, c2, trim)
                if result is not None:
                    #Only append full Dicts!
                    results.append(result)
            except ValueError:
                print(ValueError)

    result_df = pd.DataFrame(results)
    return result_df

# Parallel processing to fix!

# def pool_process(filelist, rundir, phase, ncores=2):
#     '''
#     Measures shear-wave splitting using parallel processing (Pool)
#     '''
#     runner = partial(measure_event,rundir,phase)
#     print("Parallel job. {} cores requested".format(ncores))
#     with contextlib.closing( Pool(processes = ncores) ) as pool:
#         pool.map(runner,filelist)      