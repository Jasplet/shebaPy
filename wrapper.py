#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 12:13:24 2021

@author: ja17375
"""

import os
import obspy
from obspy.taup import TauPyModel
from windower import WindowPicker
import subprocess as sub
from multiprocessing import current_process
from subprocess import CalledProcessError


class Interface:
    """
    Class which will act as the interface to sheba.
    The "subprocess" sheba will be a bound method
    """
    def __init__(self,st, phase, rundir=None):
        
        self.st = st
        self.sacstats = st[0].stats.sac
        self.station = st[0].stats.station # station code 
        self.delta = st[0].stats.delta # sample rate of seismometer [s]
        # do some basic tests of the data 
        self.check_phase_dist(phase)
        self.phase = phase # The shear-wave phase we are measuing splitting for!
        for trace in self.st:
            trace.stats.sac.kstnm = '{:>8}'.format(trace.stats.sac.kstnm)
#          Formats Station name in headers so that it is 8 characters long, with emtpy character fill with whitespaces    
            self.check_evdp(trace)
        
        self.tt_utc, self.tt_rel = self.model_traveltimes()
        if rundir is None:
            print('Setting rundir path to current working directory')
            self.path = os.getcwd()
        else:
            self.path = rundir
    
    def preprocess(self,c1=0.01,c2=0.5):
        """
        Function to bandpass filter and trim the components
        Seismograms are trimmed so that they start 1 minute before the expected arrival and end 2 minutes after the arrival
        c1 - [Hz] Lower corner frequency
        c2 - [Hz] Upper corner frequency 
        By default traces will be filtered between 0.01Hz-0.5Hz.
        Using an upper corner of 0.1Hz for SKS/SKKS is also common, but can cut out high f signal (but also reduces noise)
        """
        for trace in self.st:
#       De-mean and detrend each component
            trace.detrend(type='demean') #demeans the component
            trace.detrend(type='simple') #De-trends component
#       Filter each component. Bandpass flag gives a bandpass-butterworth filter
            trace.filter("bandpass",freqmin= c1, freqmax= c2,corners=2,zerophase=True)
#       Now trim each component to the input length
#       We only need to trim and set the window length for real data, not synthetics made with sacsplitwave
#           Now set the trim
            print('Trim Traces')
            t1 = (self.tt_utc - 60) #I.e A minute before the arrival
            t2 = (self.tt_utc + 120) #I.e Two minutes after the arrival
            trace.trim(t1,t2)
#               Initialise Windows        
            self.initialise_windows(trace)

    def measure_splitting(self,output_filename, sheba_exec_path, window=False):
        """
        Measures Shear-wave splitting using Sheba. 
        """
        fname = f'{output_filename}_{self.phase}'
        if window:
            try:
                self.window_event()
            except ValueError:
                print('Event Skipped')
                return 
            
        self.gen_infile(fname)
        self.write_out(fname)
        print(f'Worker {current_process()} Passing {fname} into Sheba.')
        out = sub.run(f'{sheba_exec_path}/sheba_exec', capture_output=True, cwd=self.rundir)

    def gen_infile(self,filename, nwind=10, tlag_max=4.0):
        '''
        Generates the input file needed for the Fortran SHEBA routines.
        '''
        with open(f'{self.path}/sheba.in','w') as writer:
            writer.write('SHEBA.IN \n')
            writer.write(f'{filename} \n') # Name of pre-processed data 
            for trace in self.st:
                writer.write('{} \n'.format(trace.stats.channel)) # Write channels for each component (E, N, Z order)
            writer.write('1 \n') # Specifies Eigenvalue minimisation, replace with spol if transverse minimisation is desired (not supported here)
            writer.write('{} {} \n'.format(nwind,nwind))
            writer.write('{} \n'.format(tlag_max)) # sets max tlag in gridsearch
            writer.write('0')  
            writer.write('0')
   
    def write_out(self,filename):
      """
      Function to write the component seismograms to SAC files to the SHEBA rundir.
      The rundir is set by self.path.
      
      filename [str] - the output filename for each trace. 
      
      """
      for trace in self.st:
          ch = trace.stats.channel
          trace.write(f'{self.path}/{filename}.{ch}', format='SAC', byteorder=1)

    def window_event(self):
        '''Function that enables manual picking of window start/end ranges'''
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
            [trace.stats.sac.update(windows) for trace in self.st]
            return
        
    def model_traveltimes(self):
        """
        Function to run TauP traveltime models for the SKS phase.
        Returns SKS predictided arrivals (seconds), origin time of the event (t0) as a UTCDateTime obejct and the SKS arrival as a UTCDateTime object
        tr - trace object for which SKS arrival time will be predicted
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

    def check_evdp(self, trace):
        '''
        Tests the event depths 
        '''
        evdp = trace.stats.sac.evdp
        if trace.stats.sac.evdp >= 1000:
            raise Warning('Event depth is greater than 1000km! EVDP may be in meters')
            trace.stats.sac({'evdp':evdp/1000})
        elif trace.stats.sac.evdp == 0:
            raise ValueError('Event depth is 0km!')

    def check_phase_dist(self,phase):
        """
        Function to test if the given phase is actually measureable!
        """
        if phase == 'SKS':
            if self.sacstats['gcarc'] > 145.0:
                raise ValueError('Event-Station distance is greater than 145, SKS not visible/reliable.')
        elif phase == 'SKKS':
            if (self.sacstats['gcarc'] < 105.0) or (self.sacstats['gcarc'] > 145.0) :
                raise ValueError('Event-Station distance outside of acceptable range of 105-145 for SKKS')
        elif phase == 'ScS':
            if self.sacstats['gcarc'] > 95.0:
                raise ValueError('Event-Station distance greatr than acceptable distance of 95 deg for ScS')
        else:
            print(f'Phase {self.phase} not ScS, SKS or SKKS. Not checking!')
