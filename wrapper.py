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
    def __init__(self,st):
        
        self.st = st
        for trace in self.st:
            trace.stats.sac.kstnm = '{:>8}'.format(trace.stats.sac.kstnm)
#          Formats Station name in headers so that it is 8 characters long, with emtpy character fill with whitespaces    
            self.check_evdp(trace)
                
        self.gcarc = st[0].stats.sac.gcarc # source-reciever great circle distance [deg]
        self.station = st[0].stats.station # station code 
        self.delta = st[0].stats.delta # sample rate of seismometer [s]
        
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
        
    def model_traveltimes(self,phase):
        """
        Function to run TauP traveltime models for the SKS phase.
        Returns SKS predictided arrivals (seconds), origin time of the event (t0) as a UTCDateTime obejct and the SKS arrival as a UTCDateTime object
        tr - trace object for which SKS arrival time will be predicted
        """
        model = TauPyModel(model="iasp91")
        tt = model.get_travel_times((self.horz1[0].stats.sac.evdp),
                                     self.horz1[0].stats.sac.gcarc,
                                     [phase])
        try:
            traveltime = tt[0].time
        except IndexError:
            print('index Error')
            traveltime = None

        evt_time = obspy.UTCDateTime(year = self.horz1[0].stats.sac.nzyear,
                                     julday = self.horz1[0].stats.sac.nzjday,
                                     hour=self.horz1[0].stats.sac.nzhour,
                                     minute=self.horz1[0].stats.sac.nzmin,
                                     second=self.horz1[0].stats.sac.nzsec,
                                     microsecond=self.horz1[0].stats.sac.nzmsec)
        self.tt = evt_time + traveltime
        self.tt_rel = traveltime
        return traveltime

    def check_phase_dist(self,phase_to_check):
        """
        Function to test if the given phase is actually measureable!
        """
        if phase_to_check == 'SKS':
            if self.gcarc > 145.0:
                raise ValueError('Event-Station distance is greater than 145, SKS not visible/reliable.')
        elif phase_to_check == 'SKKS':
            if (self.gcarc < 105.0) or (self.gcarc > 145.0) :
                raise ValueError('Event-Station distance outside of acceptable range of 105-145 for SKKS')
        elif phase_to_check == 'ScS':
            if self.gcarc > 95.0:
                raise ValueError('Event-Station distance greatr than acceptable distance of 95 deg for ScS')
        else:
            print('Phase {} not ScS, SKS or SKKS. Not checking!'.format(phase_to_check))

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

    def preprocess(self,synth=False,c1=0.01,c2=0.5):
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
            if synth == False: 
#               Now set the trim
                print('Trim Traces')
                t1 = (self.tt - 60) #I.e A minute before the arrival
                t2 = (self.tt+ 120) #I.e Two minutes after the arrival
                trace.trim(t1,t2)
#               Initialise Windows        
                self.initialise_windows(trace)

    def window_event(self):
        '''Function that enables manual picking of window start/end ranges'''
        # Windowing code
        st = self.horz1 + self.horz2
        Windower = WindowPicker(st,
                                self.st[0].stats.sac['user0'], self.st[0].stats.sac['user1'],
                                self.st[0].stats.sac['user2'], self.st[0].stats.sac['user3'],
                                self.tt_rel)
        if Windower.wbeg1 is None:
            print("Skipping")
            return False
        else:
            print("Windower Closed, adjusting window ranges")
            windows = {'user0' : Windower.wbeg1, 'user1' : Windower.wbeg2,
                      'user2' : Windower.wend1, 'user3' : Windower.wend2}
            [trace.stats.sac.update(windows) for trace in self.st]
            return True
      
    def write_out(self,phase,label,path=None):
        """
        Function to write the component seismograms to SAC files within the sheba directory structure so Sheba can access them
        station [str] - station code
        phase [str] - phase code for which seismic phase splitting is being meausured
        path [str] - path that you want the seismogrmas saved to.
        """
        for trace in self.st:
            ch = trace.stats.channel
            if path is None:
                print(f'Writing out to CWD: {path}')
                path = os.get_cwd()
                
            trace.write(f'{path}/{label}{phase}.{ch}', format='SAC', byteorder=1)

    def gen_infile(self,path,label,phase,nwind=10,tlag_max=4.0):
    
        os.chdir(path) # Make sure we are in the right directory
    
        with open('sheba.in','w') as writer:
            writer.write('SHEBA.IN \n')
            writer.write('{}{} \n'.format(label,phase)) # write file prefix
            for trace in self.st:
                writer.write('{} \n'.format(trace.stats.channel)) # Write channels for each component (E, N, Z order)
            writer.write('1 \n') # Specifies Eigenvalue minimisation, replace with spol if transverse minimisation is desired (not supported here)
            writer.write('{:i} {:i} \n'.format(nwind,nwind))
            writer.write('{} \n'.format(tlag_max)) # sets max tlag in gridsearch
            writer.write('0')  
            writer.write('0')


    def sheba(self,phase,label,path=None,nwind=True,):
        """
        The big one! This function uses the subprocess module to host sac and then runs sheba as a SAC macro
        """
        print(f'Worker {current_process()} Passing {label} into Sheba for.')
        p = sub.Popen(['sac'],
                     stdout = sub.PIPE,
                     stdin  = sub.PIPE,
                     stderr = sub.STDOUT,
                     cwd = path,
                     encoding='utf8')
#       Now that the child process SAC has been opened, lets specify what arguements we want to pass to it
#       echo on means commands should be echoed (so we can check whats going on if we print output)
#       SETMACRO makes sure SAC can find sheba (avoids pathing problems)
#       m sheba calls sheba as a SAC macro
        if nwind == True:
            s = '''
            echo on\n
            SETMACRO /Users/ja17375/Ext_programs/macros
            m sheba file {}{} plot yes pick no nwind 10 10 batch yes
            '''.format(label,phase)
        elif nwind == False:
            s = '''
            echo on\n
            SETMACRO /Users/ja17375/Ext_programs/macros
            m sheba file {}{} plot yes pick yes batch yes
            '''.format(label,phase)
        try:           
            _ = p.communicate(s)
            # print(out[0])
        except CalledProcessError as err:
            print(err)