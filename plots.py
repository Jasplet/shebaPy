#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 14:45:47 2021

@author: ja17375

Functions to plot shear-wave ata and shear-wave splitting results
"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import obspy
from obspy import UTCDateTime
import netCDF4 as nc

def diagnostic_plot(result_nc, result_id):
    '''
    Produces a diagnostic plot for each shear-wave splitting measurment

    Shows:
        - corrected/uncorrected traces, with windows and S pick
        - corrected/uncorrected particle motions
        - normalised egigenvalue surface
    '''


    fig = plt.figure(layout="constrained")

    gs = GridSpec(3, 4, figure=fig)

    # Add uncorrected traces
    ax1 = fig.add_subplot(gs[0,0:2])

    ax2 = fig.add_subplot(gs[1,0:2])
    ax3 = fig.add_subplot(gs[2, 0])
    ax4 = fig.add_subplot(gs[2, 1])
    ax4 = fig.add_subplot(gs[1:3, 2:4])

def _plot_traces(st, **kwargs):
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
    
    if 'show_final_window' in kwargs:
        ax.axvline(x=st[0].stats.sac['a'], linewidth=1, color='black', linestyle='-')
        ax.axvline(x=st[0].stats.sac['f'], linewidth=1, color='black', linestyle='-')
    
    # set axis label
    if 'units' not in kwargs:
        kwargs['units'] = 's'     
            
    ax.set_xlabel(f'Time relative to origin ({kwargs["units"]})')
    ax.legend(framealpha=0.75)
    
    return 

def _ppm(ax, st, wbeg, wend):
    st_plot = st.copy()
    
    rel_msecs = int(st[0].stats.sac['nzmsec']*1e-3)
    rel_secs = int(st[0].stats.sac['nzsec'])
    rel_mins = int(st[0].stats.sac['nzmin'])
    rel_hours= int(st[0].stats.sac['nzhour'])
    rel_jdays = int(st[0].stats.sac['nzjday'])
    rel_years = int(st[0].stats.sac['nzyear'])
    event_time = obspy.UTCDateTime(year=rel_years, julday=rel_jdays, hour=rel_hours,
                                            minute=rel_mins, second=rel_secs, microsecond=rel_msecs)
    st_plot.trim(event_time + wbeg, event_time + wend)
    trN = st_plot.select(channel='??N')[0]
    trE = st_plot.select(channel='??E')[0]
    ax.plot(trE.data, trN.data)
    ax.set_xlabel('East')
    ax.set_ylabel('North')
    return 


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