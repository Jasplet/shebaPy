#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 12:13:24 2021

@author: ja17375
"""

import os
import subprocess as sub
from netCDF4 import Dataset
import obspy
from obspy.taup import TauPyModel
from obspy import UTCDateTime
from windower import WindowPicker

# from .plots import plot_traces, plot_pm

class Wrapper:
    """
    Class that 'wraps' around SHEBA fortran routines.
    
    Includes methods for pre-processing and windowing shear-wave data
    
    Attributes
    ----------
    st :
        obspy Stream object that holds waveform data
    sacstats : dict
        metadata from SAC headers
    station : str 
        Recording station code
    delta : float
        seismometer sample rate
    phase : str
        seismic phase recorded
    tt_utc :
        traveltime in UTCDateTime format
    tt_rel :
        relative traveltime (i.e., seconds after event time)
    path : str
        path to working directory for SHEBA
    Methods
    ----------
    
    
    """
    def __init__(self,st, phase, rundir=None, **kwargs):
        '''
        Constructs Wrapper for a 3-component waveform
        
        Parameters:
        ----------
        st : 
            obspy Stream object containing raw data
        phase : str
            seismic phase of interes
        rundir : str, optional, default=None
            path to run directory. if None current working directory is used
        '''
        self.st = st
        self.sacstats = st[0].stats.sac
        self.station = st[0].stats.station # station code
        self.delta = st[0].stats.delta # sample rate of seismometer [s]
        self.phase = phase # The shear-wave phase we are measuing splitting for!
        for trace in self.st:
            trace.stats.sac.kstnm = '{:>8}'.format(trace.stats.sac.kstnm)
            self.fix_cmp_dir(trace)
#       Formats Station name in headers so that it is 8 characters long,
#       with emtpy character fill with whitespaces

        if 'wbeg_pre_S' in kwargs:
            self.wbeg_pre_S = kwargs['wbeg_pre_S']
        else: 
            self.wbeg_pre_S = 0.1

        if 'wend_post_S' in kwargs:
            self.wend_post_S = kwargs['wend_post_S']
        else:
            self.wend_post_S = 0.2
        if 'pick_tol' in kwargs:
            self.pick_tol = kwargs['pick_tol']
        else:
            self.pick_tol = 0.05

        check_phase_dist(phase, self.sacstats['gcarc'])
        check_evdp(trace)
        self.tt_utc, self.tt_rel = self.model_traveltimes()
        if rundir is None:
            print('Setting rundir path to current working directory')
            self.path = os.getcwd()
        else:
            self.path = rundir

    def fix_cmp_dir(self, trace):
        """
        Fixes the sac headers cmpinc, cmpaz for BHE,BHN, BHZ channels. If data is 
        downlaoded directly from the IRIS DMC then these SAC headers are missing
        """
        if trace.stats.channel == 'BHE':
            trace.stats.sac.cmpinc = 90
            trace.stats.sac.cmpaz = 90
        elif trace.stats.channel == 'BHN':
            trace.stats.sac.cmpinc = 90
            trace.stats.sac.cmpaz = 0
        elif trace.stats.channel == 'BHZ': 
            trace.stats.sac.cmpinc = 0
            trace.stats.sac.cmpaz = 0

    def preprocess(self,c1=0.01,c2=0.5, trim=True):
        """
        Function to bandpass filter and trim the components
        Seismograms are trimmed so that they start 1 minute before the expected arrival 
        and end 2 minutes after the arrival.
        By default traces will be filtered between 0.01Hz-0.5Hz.
        Using an upper corner of 0.1Hz for SKS/SKKS is also common, 
        but can cut out high f signal (but also reduces noise)
        
        Parameters
        ----------
        c1 : float, optional, default=0.01
            Lower corner frequency [Hz]
        c2 : float, optional, default=0.5
            Upper corner frequency [Hz]
 
        """
        for trace in self.st:
#       De-mean and detrend each component
            trace.detrend(type='demean') #demeans the component
            trace.detrend(type='simple') #De-trends component
#       Filter each component. Bandpass flag gives a bandpass-butterworth filter
            trace.filter("bandpass",freqmin= c1, freqmax= c2,corners=2,zerophase=True)
#       Data is only trimmed if the record is longer than 3 minutes.
#       This is to ensure there is enough space for windowing, especially manual windowing. 
#       This record length makes sense for teleseismic SKS, SKKS, ScS data
#       but may need revising for other shear-wave data.
            if ((trace.stats.endtime - trace.stats.starttime) > 180.0) & trim == True:
                t1 = (self.tt_utc - 60) #I.e A minute before the arrival
                t2 = (self.tt_utc + 120) #I.e Two minutes after the arrival
                trace.trim(t1,t2)
            else:
                print('Not trimming')
#               Initialise Windows
            self.initialise_windows(trace)

    def measure_splitting(self,output_filename, sheba_exec_path, window=False, nwind=10, debug=False, **kwargs):
        """
        Measures Shear-wave splitting using Sheba. 

        Parameters
        ----------
        output_filename : str
            filestem to name processed data and SHEBA outfiles
        sheba_exec_path : str
            path to the compiled SHEBA executable sheba_exec
        window : bool, optional, default=False
            switch to manual windowing (if True) or to use pre-defined windows (False)
        nwind : int, optional, default=10
            number of window start/ends to consider in SHEBA cluster analysis. nwind=10 tries 100 combinations
        debug : bool, optional, default=False
            prints out SHEBA stdout when True. Default is False
            
        Returns:
            result : dict
                Shear-wave splitting measurement results and metadata
        """
        if window:
            try:
                #print(self.st)
                self.window_event()
            except ValueError:
                print('Event Skipped')
                return
        else:
            self.auto_window()
        if 'tlag_max' in kwargs:
            self.gen_infile(output_filename, nwind=nwind, tlag_max=kwargs['tlag_max'])
        else:
            self.gen_infile(output_filename, nwind=nwind)

        self.write_out(output_filename)
        #print(f'Passing {output_filename} into Sheba.')
        out = sub.run(f'{sheba_exec_path}/sheba_exec', capture_output=True, cwd=self.path, check=True)
        if debug:
            # print what sheba returns to stdout. useful for debugging the wrapping.
            print(out)
        result = collate_result(self.path, output_filename)
        self.update_sachdrs(output_filename, result)
        return result
        
    def update_sachdrs(self, filename, result):
        '''
        Updates the sac headers of traces to add measured windows 
        Parameters
        ----------
        result : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        for trace in self.st:
            trace.stats.sac.update({'a':result['WBEG'], 'f':result['WEND']})
        self.write_out(filename)
        #Also 
        st_corr = obspy.read(f'{self.path}/{filename}_corr.?H?')
        for trace in st_corr:
            ch = trace.stats.channel
            trace.stats.sac.update({'a':result['WBEG'], 'f':result['WEND']})
            trace.write(f'{self.path}/{filename}_corr.{ch}', format='SAC', byteorder=1)
                
    def gen_infile(self,filename, nwind=10, tlag_max=4.0):
        '''
        Generates the input file needed for the Fortran SHEBA routines.
        
        Parameters
        ----------
        filename : str
            name of waveform data for sheba to read in 
        nwind : int, optional, default=10
            number of analysis window starts/ends to use in cluster analysis
        tlag_max : flat, optional, default=4.0
            maximum allowable delay time in grid search
            
        '''
        with open(f'{self.path}/sheba.in','w') as writer:
            writer.write('SHEBA.IN \n')
            writer.write(f'{filename} \n') # Name of pre-processed data
            for trace in self.st:
                writer.write('{} \n'.format(trace.stats.channel))
            # 1 Specifies Eigenvalue minimisation, 
            # replace with spol if transverse minimisation is desired (not supported here)
            writer.write('1 \n')
            writer.write(f'{nwind} {nwind} \n')
            writer.write(f'{tlag_max} \n') # sets max tlag in gridsearch
            writer.write('0 \n')
            writer.write('0')

    def write_out(self,filename):
        """
        Function to write the component seismograms to SAC files to the SHEBA rundir.
        The rundir is set by self.path.
        
        Parameters
        -----------
        filename : str
            the output filename for each trace
        
        """
        for trace in self.st:
            ch = trace.stats.channel
            trace.write(f'{self.path}/{filename}.{ch}', format='SAC', byteorder=1)

    def window_event(self, **kwargs):
        '''
        Function that enables manual picking of window start/end ranges by initialising a 
        WindowPicker object
        '''
        # Windowing code
        Windower = WindowPicker(self.st,
                                self.st[0].stats.sac['user0'], self.st[0].stats.sac['user1'],
                                self.st[0].stats.sac['user2'], self.st[0].stats.sac['user3'],
                                self.tt_rel)

        if Windower.wbeg1 is None:
            raise ValueError('Skipping event as it is poor quality!')
        else:
            print("Windower Closed, adjusting window ranges")
            windows = {'user0' : Windower.wbeg1, 'user1' : Windower.wbeg2,
                      'user2' : Windower.wend1, 'user3' : Windower.wend2}
            for trace in self.st:
                trace.stats.sac.update(windows)
            return

    def model_traveltimes(self):
        """
        Uses TauP to predict phase traveltime.
        
        Returns
        ----------
        tt_utc :
            predicted arrival time as a UTCDateTime object
        traveltime :
            predicted arrival time as seconds after event origin time
        """
        evt_time = obspy.UTCDateTime(year = self.sacstats['nzyear'],
                                     julday = self.sacstats['nzjday'],
                                     hour=self.sacstats['nzhour'],
                                     minute=self.sacstats['nzmin'],
                                     second=self.sacstats['nzsec'],
                                     microsecond=self.sacstats['nzmsec'])
        if ('t1' in self.sacstats) & (self.phase == 'SKS'):
            # pick exists and is stored in sac headers
            print('Use existing pick in t1')
            traveltime = self.sacstats['t1']
        else:
            print(f'Use TauP to predict {self.phase} arrival time') 
            model = TauPyModel(model="iasp91")
            tt = model.get_travel_times((self.sacstats['evdp']),
                                         self.sacstats['gcarc'],
                                         [self.phase])  
            traveltime = tt[0].time
        
        tt_utc =  evt_time + traveltime + self.st[0].stats.sac['o']
        print(self.st[0].stats.sac['o'])

        print(f'Depth: {self.sacstats["evdp"]}, Epicentral distance: {self.sacstats["gcarc"]}')
        print(tt_utc)
        print(traveltime)
        return tt_utc, traveltime

    def initialise_windows(self, trace):
        '''
        Checks for the presence of analysis window in trace SAC headers (user0-3).
        If they are not present default values are initialised around the predicted arrival time
        
        Parameters:
        ----------
        trace : 
            obspy Trace object to test
        '''
        if all (k in trace.stats.sac for k in ('user0','user1','user2','user3')):
            user0 = trace.stats.sac['user0']
            user1 = trace.stats.sac['user1']
            user2 = trace.stats.sac['user2']
            user3 = trace.stats.sac['user3']
        else:
            if 't2' in trace.stats.sac:
                user0, user1, user2, user3 = auto_window(trace.stats.sac['t2'], self.wbeg_pre_S, self.wend_post_S, pick_tol=self.pick_tol)
            else:
                user0, user1, user2, user3 = auto_window(self.tt_rel, wbeg_pre_S=15, wend_post_S=15)    

            keychain = {'user0':user0,'user1':user1,'user2':user2,'user3':user3}
            trace.stats.sac.update(keychain)

def collate_result(path=None, fname=None, full_file=None):
        '''
        Collates meausurement results from SHEBA, allowing them to be used by other functons
        
        Parameters
        ----------  
        path : str
            path to run directory with NetCDFs result files to read.

        fname : str
            output filestem used by SHEBA
        
        Returns
        ----------
        result : dict 
            Shear-wave splitting measurement results and metadata
        '''
        if full_file is None:
            raw_result = Dataset(f'{path}/{fname}_sheba_result.nc')
        elif (path is None) and (fname is None):
            print(full_file)
            raw_result = Dataset(full_file)
        # print('Best fitting result is')
        # print(f'Fast direction =  {raw_result.fast} +/- {raw_result.dfast}')
        # print(f'Delay time = {raw_result.tlag} +/- {raw_result.dtlag}')

        result = {'STAT':raw_result.station.strip(),
                'DATE':raw_result.zdate,'TIME':raw_result.ztime.split('.')[0],
                'STLA':raw_result.stla, 'STLO':raw_result.stlo,
                'EVLA':raw_result.evla, 'EVLO':raw_result.evlo, 'EVDP':raw_result.evdp,
                'GCARC':raw_result.gcarc, 'AZI':raw_result.az, 
                'BAZ':raw_result.baz, 'SPOL':raw_result.spol,
                'WBEG':raw_result.wbeg, 'WEND':raw_result.wend, 
                'FAST':raw_result.fast, 'DFAST':raw_result.dfast,
                'TLAG':raw_result.tlag, 'DTLAG':raw_result.dtlag,
                'SI(Pa)':raw_result.intensity_estimated,
                'SI(Pr)':raw_result.intensity, 'Q':raw_result.qfactor,
                'EIGORIG':raw_result.eigrat_orig, 'EIGCORR': raw_result.eigrat_corr,
                'SNR':raw_result.snr, 'NDF':raw_result.ndf}
        return result

def event_relative_time(st_stats):
    '''
    Use SAC headers
    
    '''
    start = st_stats.starttime 
    sacstat = st_stats.sac
    
    startdate = UTCDateTime(start)
    eventtime = UTCDateTime(year=sacstat['nzyear'], julday=sacstat['nzjday'], hour=sacstat['nzhour'],
                            minute=sacstat['nzmin'], second=sacstat['nzsec'], microsecond=sacstat['nzmsec'])
    rel_start = startdate - eventtime
    return rel_start, eventtime

def check_evdp(trace):
    '''
    Tests the event depth in the trace SAC headers to make sure it is a sensible value
    
    Parameters:
    ----------
    trace : 
        obspy Trace object to test    
    '''
    evdp = trace.stats.sac.evdp
    if trace.stats.sac.evdp >= 1000:
        trace.stats.sac({'evdp':evdp/1000})
        raise Warning('Event depth is greater than 1000km! EVDP may be in meters')
    elif trace.stats.sac.evdp == 0:
        raise ValueError('Event depth is 0km!')

def check_phase_dist(phase, gcarc):
    """
    Function to test if the given phase is actually measureable
    
    Parameters
    ----------
    phase : str 
        seismic phase to measure
    gcarc : float
        great circle distance between earthquake and the recording station
    """
    if phase == 'SKS':
        if gcarc > 145.0:
            raise ValueError(f'Event-Station distance {gcarc} is greater than 145, SKS not visible/reliable.')
    elif phase == 'SKKS':
        if (gcarc < 105.0) or (gcarc > 145.0) :
            raise ValueError(f'Event-Station distance {gcarc} outside of acceptable range of 105-145 for SKKS')
    elif phase == 'ScS':
        if gcarc > 95.0:
            raise ValueError(f'Event-Station distance {gcarc} greater than acceptable distance of 95 deg for ScS')
    elif phase == 'Synth':
        print('Skip test for synthetics')
    else:
        print(f'Phase {phase} not ScS, SKS or SKKS. Not checking!')

def auto_window(S_pick, wbeg_pre_S=0.1, wend_post_S=0.2, pick_tol=0.05):
        '''
        Automatically window traces based on an assumed S pick
        
        Assumes S pick is in t2 in SAC header
        '''
        wbeg1 = S_pick - pick_tol - wbeg_pre_S
        wbeg2 = S_pick - pick_tol
        wend1 = S_pick + pick_tol + wend_post_S
        wend2 = S_pick + pick_tol + wend_post_S + (wbeg2 - wbeg1)

        return wbeg1, wbeg2, wend1, wend2