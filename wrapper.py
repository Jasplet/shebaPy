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
    def __init__(self,st, phase, rundir=None):
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
#       Do some basic data QA
        if phase != 'Synth':
            print(phase)
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

    def preprocess(self,c1=0.01,c2=0.5):
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
            if (trace.stats.endtime - trace.stats.starttime) > 180.0:
                print('Trim Traces')
                t1 = (self.tt_utc - 60) #I.e A minute before the arrival
                t2 = (self.tt_utc + 120) #I.e Two minutes after the arrival
                trace.trim(t1,t2)
            else:
                print('Traces too short (<180 s ) to trim')
#               Initialise Windows
            self.initialise_windows(trace)

    def measure_splitting(self,output_filename, sheba_exec_path, window=False, debug=False):
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
        debug : bool, optional, default=False
            prints out SHEBA stdout when True. Default is False
            
        Returns:
            result : dict
                Shear-wave splitting measurement results and metadata
        """
        if window:
            try:
                self.window_event()
            except ValueError:
                print('Event Skipped')
                return

        self.gen_infile(output_filename)
        self.write_out(output_filename)
        print(f'Passing {output_filename} into Sheba.')
        out = sub.run(f'{sheba_exec_path}/sheba_exec', capture_output=True, cwd=self.path, check=True)
        if debug:
            # print what sheba returns to stdout. useful for debugging the wrapping.
            print(out)
        result = self.collate_result(output_filename)
        return result

    def collate_result(self, fname):
        '''
        Collates meausurement results from SHEBA, allowing them to be used by other functons
        
        Parameters
        ----------    
        fname : str
            output filestem used by SHEBA
        
        Returns
        ----------
        result : dict 
            Shear-wave splitting measurement results and metadata
        '''

        raw_result = Dataset(f'{self.path}/{fname}_sheba_result.nc')
        print('Best fitting result is')
        print(f'Fast direction =  {raw_result.fast} +/- {raw_result.dfast}')
        print(f'Delay time = {raw_result.tlag} +/- {raw_result.dtlag}')

        result = {'STAT':raw_result.station.strip(),
                'DATE':raw_result.zdate,'TIME':raw_result.ztime.split('.')[0],
                'STLA':raw_result.stla, 'STLO':raw_result.stlo,
                'EVLA':raw_result.evla, 'EVLO':raw_result.evlo, 'EVDP':raw_result.evdp,
                'GCARC':raw_result.gcarc, 'AZI':raw_result.az, 'BAZ':raw_result.baz, 
                'WBEG':raw_result.wbeg, 'WEND':raw_result.wend, 
                'FAST':raw_result.fast, 'DFAST':raw_result.dfast,
                'TLAG':raw_result.tlag, 'DTLAG':raw_result.dtlag,
                'SI(Pa)':raw_result.intensity_estimated,
                'SI(Pr)':raw_result.intensity, 'Q':raw_result.qfactor,
                'SNR':raw_result.snr, 'NDF':raw_result.ndf}
        return result
        
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

    def window_event(self):
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
        model = TauPyModel(model="iasp91")
        tt = model.get_travel_times((self.sacstats['evdp']),
                                     self.sacstats['gcarc'],
                                     [self.phase])  
        traveltime = tt[0].time
        evt_time = obspy.UTCDateTime(year = self.sacstats['nzyear'],
                                     julday = self.sacstats['nzjday'],
                                     hour=self.sacstats['nzhour'],
                                     minute=self.sacstats['nzmin'],
                                     second=self.sacstats['nzsec'],
                                     microsecond=self.sacstats['nzmsec'])
        tt_utc =  evt_time + traveltime
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
            print('Pre-defined window ranges found')
        else:
            user0 = self.tt_rel - 15 # 15 seconds before arrival
            user1 = self.tt_rel # t predicted arrival
#           Set the range of window endtime (user2/user3)
            user2 = self.tt_rel + 15 # 15 seconds after, gives a min window size of 20 seconds
            user3 = self.tt_rel + 30 # 30 seconds after, gives a max window size of 45 seconds
            keychain = {'user0':user0,'user1':user1,'user2':user2,'user3':user3}
            trace.stats.sac.update(keychain)

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
