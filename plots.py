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

def diagnostic_plot(st, st_corr, result, event_time):
    '''
    Produces a diagnostic plot for each shear-wave splitting measurment

    Shows:
        - corrected/uncorrected traces, with windows and S pick
        - corrected/uncorrected particle motions
        - normalised egigenvalue surface
    '''


    fig = plt.figure(layout="constrained", figsize= (13,7.5))
    gs = GridSpec(4, 6, figure=fig)
    # Input data
    ax1 = fig.add_subplot(gs[0,0:3])
    _plot_traces(st, show_final_window=True, axes=ax1)
    ax1.set_xlim([result.wbeg -1, result.wend + 1])
    ax1.set_title(f'Input S. Event: {st[0].stats.starttime}. Station: {result.station.strip()}', fontsize=12)
    # Corrected data
    ax2 = fig.add_subplot(gs[0,3:], sharey=ax1)
    _plot_traces(st_corr, show_final_window=True, axes=ax2)
    ax2.set_xlim([result.wbeg -1, result.wend + 1])
    ax2.set_title(f'Corrected S. $\phi_f = {result.fast:4.2f}\pm{result.dfast:4.2f}$, $\delta t = {result.tlag:4.3f}\pm{result.dtlag:4.3f}$')
    # 
    ax3 = fig.add_subplot(gs[1, 0])
    _ppm(ax3, st, event_time)
    
    ax4 = fig.add_subplot(gs[1, 1])
    _ppm(ax4, st_corr, event_time)
    ax4 = fig.add_subplot(gs[1:4, 3:6])
    phis = result.variables['fast_vector'][:]
    dts = result.variables['tlag_vector'][:]
    PHI, TLAG = np.meshgrid(phis, dts)

    C = ax4.contourf(TLAG, PHI, result.variables['lam2_norm_grid'][:], levels=np.linspace(1,21,11),
                cmap='magma_r')
    C1 = ax4.contour(TLAG, PHI, result.variables['lam2_norm_grid'][:], levels=[1], linewidths=[3],
                colors='black')
    ax4.plot(result.tlag, result.fast, 'x', color='royalblue')
    ax4.set_ylabel('Fast polarisation [°]')
    ax4.set_xlabel('Delay time [s]')
    ax4.set_xlim([dts.min(), dts.max()])
    ax4.set_ylim([phis.min(), phis.max()])

    # Window clustser fast
    wind_num = np.arange(0,result.dimensions['window'].size) + 1 #shift to start from 1
    ax5 = fig.add_subplot(gs[2,0:3])
    ax5.errorbar(x=wind_num, y=result.variables['mw_fast'][:], yerr=result.variables['mw_dfast'][:])
    ax5.errorbar(x=wind_num[result.best_window-1], y=result.variables['mw_fast'][result.best_window-1],
                yerr=result.variables['mw_dfast'][result.best_window-1], color='tab:red')
    ax5.set_ylim([phis.min(),phis.max()])
    ax5.set_ylabel(r'$\phi_f$ [°]')
    # dt for all windows
    ax6 = fig.add_subplot(gs[3,0:3])
    ax6.errorbar(x=wind_num, y=result.variables['mw_tlag'][:], yerr=result.variables['mw_dtlag'][:])
    ax6.errorbar(wind_num[result.best_window-1], y=result.variables['mw_tlag'][result.best_window-1],
              yerr=result.variables['mw_dtlag'][result.best_window-1], color='tab:red')
    ax6.set_ylim([dts.min(), dts.max()])
    ax6.set_ylabel(r'$\delta t$ [s]')
    ax6.set_xlabel('Window #')
    # Clusters
    ax7 = fig.add_subplot(gs[1,2])
    ax7.scatter(x=result.variables['cluster_xc0'][:], y=result.variables['cluster_yc0'][:])
    ax7.plot(result.variables['cluster_xc0'][result.best_cluster-1], result.variables['cluster_yc0'][result.best_cluster-1],'x', color='red')
    ax7.set_xlim([dts.min(), dts.max()])
    ax7.set_ylim([phis.min(), phis.max()])
    ax7.set_ylabel(r'$\delta t$ [s]')
    ax7.set_ylabel(r'$\phi_f$ [°]')

    return fig

def _plot_traces(st, event_time, **kwargs):
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
    
    times = st.times(reftime=event_time)
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

def _ppm(ax, st, event_time):
    st_plot = st.copy()

    st_plot.trim(event_time + st[0].stats.sac['a'], event_time + st[0].stats.sac['f'])
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
    if 'o' in st[0].stats.sac:
        relative_times = st[0].times() + st[0].stats.sac['b'] + st[0].stats.sac['o']
    else:
        relative_times = st[0].times() + st[0].stats.sac['b'] 
    return relative_times
    
if __name__ == '__main__':
    st = obspy.read('/Users/eart0593/Projects/SHARP/splitting/data/ML_gte_0/*/GB_MONM_20180420145936.*')

    _plot_traces(st)