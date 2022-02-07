#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 14:45:47 2021

@author: ja17375

Functions to plot shear-wave ata and shear-wave splitting results
"""

import matplotlib.pyplot as plt
import numpy as np
import obspy
from obspy import UTCDateTime

def plot_traces(st, **kwargs):
    '''
    function to plot shear-wave traces 
    
    Parameters:
    ----------
    st :
        obspy Stream conatining waveform data to plot
    '''
    
    if 'axes' not in kwargs:
        # no axes provided so make our own
        fig, ax = plt.subplots(1, 1)
    else:
        ax = kwargs['axes']
    
    times = time_shift(st)
    
    ax.plot(times, st[0].data, label=st[0].stats.channel, color='dodgerblue')
    ax.plot(times, st[1].data, label=st[1].stats.channel, color='darkorange')
    
    if 'show_window_range' in kwargs:
        # window range should be a list of [wbeg1 wend1 wbeg2 wend2]
        for marker in ['user0', 'user1', 'user2', 'user3']:
            ax.axvline(x=st[0].stats.sac[marker], linewidth=1, color='black', linestyle='--')
        
    
    # set axis label
    if 'units' not in kwargs:
        kwargs['units'] = 's'     
            
    ax.set_xlabel(f'Time relative to origin ({kwargs["units"]})')
    ax.legend(framealpha=0.75)
    plt.show()

def time_shift(st):
    '''
    Calculates time shift to make st.times() in seconds after event 
    
    Parameters:
    ----------
    st :
        obspy Stream conatining waveform data
        
    Returns
    ----------
    relative_times : array-like
        event-relative times (i.e., time after earthquake source time)
    '''
    sacstats = st[0].stats.sac
    origin = UTCDateTime(year=sacstats['nzyear'], julday=sacstats['nzjday'],
                               hour=sacstats['nzhour'], minute=sacstats['nzmin'],
                               second=sacstats['nzsec'], microsecond=sacstats['nzmsec']
                              )
    time_after_origin = st[0].stats.starttime - origin
    relative_times = st[0].times() + time_after_origin
    return relative_times
    
if __name__ == '__main__':
    st = obspy.read('/Users/ja17375/Programs/shebaPy/example/HKT_example*.BH?')
    plot_traces(st)