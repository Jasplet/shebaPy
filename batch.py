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
SHEBA_EXEC = '/Users/ja17375/Ext_Programs/bin'

def measure_event(file, rundir, phase, window=False, plot=False, nwind=10, debug=False, c1=0.01, c2=0.5, trim=True, **kwargs):
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
    if window:
        print('Manual windowing enabled')
    else:
        print('Will set auto window params')
        if 'wbeg_pre_S' in kwargs:
            wbeg_pre_S = kwargs['wbeg_pre_S']
        else: 
            wbeg_pre_S = 0.1
        if 'wend_post_S' in kwargs:
            wend_post_S = kwargs['wend_post_S']
        else:
            wend_post_S = 0.2
        if 'pick_tol' in kwargs:
            pick_tol = kwargs['pick_tol']
        else:
            pick_tol = 0.05 
        print(f'Earliest window start: {wbeg_pre_S:4.2f}s pre S pick')
        print(f'Earliest window end: {wend_post_S:4.2f}s post S pick')
        print(f'Pick tolerence {pick_tol:4.2f}s')

    input_data = obspy.read(file)
    sheba_wrapper = Wrapper(input_data, phase, rundir=rundir)
    sheba_wrapper.preprocess(c1, c2, trim)
    outfile = file.split('.')[0].split('/')[-1]
    # selects just filename with not extentsion or preceding path
    if 'tlag_max' in kwargs:
        splitting_result = sheba_wrapper.measure_splitting(outfile, SHEBA_EXEC, window, nwind, debug, tlag_max=kwargs['tlag_max'])
    else:
        splitting_result = sheba_wrapper.measure_splitting(outfile, SHEBA_EXEC, window, nwind, debug)
    if plot:
        print('Plotting to be added soon!')
    return splitting_result

def serial_process(filelist, rundir, phases, window=False, nwind=10, debug=False, **kwargs):
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
    #some defaults
    c1=0.01
    c2=0.5
    trim=True
    remeasure=False, 
    tlag_max = 4.0
    if 'c1' in kwargs:
        c1 = kwargs['c1']
    if 'c2' in kwargs:
        c2 = kwargs['c2']
    if 'trim' in kwargs:
        trim = kwargs['trim']
    if 'remeasure' in kwargs:
        remeasure = kwargs['remeasure']
    if 'tlag_max' in kwargs:
        tlag_max = kwargs['tlag_max']

    results = []


    for i, file in enumerate(filelist):
        res_stem = file.split('/')[-1]
        res_id = res_stem.split('.')[0]
        if os.path.isfile(f'{rundir}/{res_id}_sheba_result.nc'):
            if remeasure:
                try:
                    result = measure_event(file, rundir, phases[i], window, nwind, debug, c1, c2, trim, tlag_max=tlag_max)
                    if result is not None:
                        #Only append full Dicts!
                        results.append(result)
                except ValueError:
                    print(ValueError)
            else: 
                print('Skipping event - already measured')
                result = collate_result(rundir, res_id)
        else:
            try:
                result = measure_event(file, rundir, phases[i], window, nwind, debug, c1, c2, trim, tlag_max=tlag_max)
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