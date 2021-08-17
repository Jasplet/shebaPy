#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 17 12:13:24 2021

@author: ja17375
"""
class Interface:
    """
    Class which will act as the interface to sheba.
    The "subprocess" sheba will be a bound method
    """
    def __init__(self,st):
        # self.date = date
        # self.time = time
        self.ch = 'H'
        for i in [0,1,2]:
            st[i].stats.sac.kstnm = '{:>8}'.format(st[i].stats.sac.kstnm)
#          Formats Station name in headers so that it is 8 characters long, with emtpy character fill with whitespaces
        try:
            self.BHE = st.select(channel='BHE')
            self.BHE[0].stats.sac.cmpinc = 90
            self.BHE[0].stats.sac.cmpaz = 90
            self.BHN = st.select(channel='BHN')
            self.BHN[0].stats.sac.cmpinc = 90
            self.BHN[0].stats.sac.cmpaz = 0
            self.BHZ = st.select(channel='BHZ')
            self.BHZ[0].stats.sac.cmpinc = 0
            self.BHZ[0].stats.sac.cmpaz = 0

        except IndexError:
            try:
                print('BH? channels not found, trying BX?')
                self.BHE = st.select(channel='BXE')
                self.BHE[0].stats.sac.cmpinc = 90
                self.BHE[0].stats.sac.cmpaz = 90
                self.BHN = st.select(channel='BXN')
                self.BHN[0].stats.sac.cmpinc = 90
                self.BHN[0].stats.sac.cmpaz = 0
                self.BHZ = st.select(channel='BXZ')
                self.BHZ[0].stats.sac.cmpinc = 0
                self.BHZ[0].stats.sac.cmpaz = 0
                self.ch = 'X'
            except IndexError:
                    print('Channel selection failed, indexing instead')
                    self.BHE = st[0]
                    self.BHE.stats.sac.cmpinc = 90
                    self.BHE.stats.sac.cmpaz = 90
                    self.BHN = st[1]
                    self.BHN.stats.sac.cmpinc = 90
                    self.BHN.stats.sac.cmpaz = 0
                    self.BHZ = st[2]
                    self.BHZ.stats.sac.cmpinc = 0
                    self.BHZ.stats.sac.cmpaz = 0
        # print(self.BHN)
#       Also lets load the gcarc from each stream, so we can test for whether SKKS should be measuable
        print('Checking for lost keys')
        if all (k in self.BHE[0].stats.sac for k in ('user0','user1','user2','user3')):
            print('No lost keys!')
        else:
            print('We\'ve lost at least one key, added to dict')
            keychain = {'user0':0,'user1':1,'user2':2,'user3':3}
            self.BHE[0].stats.sac.update(keychain)
            self.BHN[0].stats.sac.update(keychain)
            self.BHZ[0].stats.sac.update(keychain)
        self.gcarc = (st[0].stats.sac.gcarc)
        self.station = st[0].stats.station
        self.delta = st[0].stats.delta
        self.bad = False
        # self.path = os.getcwd()
#       As this gcarc is calculated in split_read.py I know that it should be the same for all three traces
#       So for ease we will always read it from st[0]

    def model_traveltimes(self,phase):
        """
        Function to run TauP traveltime models for the SKS phase.
        Returns SKS predictided arrivals (seconds), origin time of the event (t0) as a UTCDateTime obejct and the SKS arrival as a UTCDateTime object
        tr - trace object for which SKS arrival time will be predicted
        """
        model = ob.taup.tau.TauPyModel(model="iasp91")

        # Add in a test for the case where depth has been gien in meters (as OBSPY IS DUMB AS HELL AND RETURNS DEPTH IN [m] FFS)
        if self.BHN[0].stats.sac.evdp >= 1000:
            #This might be a bit generous but my events shouldnt be shallower the ~ 1 km or deeper than 1000km anyway (and if this is the case than there is something SERIOUSLY wrong with our earth models)
            traveltime = model.get_travel_times((self.BHN[0].stats.sac.evdp/1000),self.BHN[0].stats.sac.gcarc,[phase])[0].time
            print(traveltime)
        elif self.BHN[0].stats.sac.evdp == 0: # Theres an event where the event data couldnt be found so evdp was set to be 0
            # Having a depth of zero will give us problems so NOW change it to 10.0km exactly (these traveltimes could be very dodgy)
            err_out = open('/Users/ja17375/DiscrePy/Sheba/Events_with_evdp_of_0.txt','w+')
            err_out.write('Station: {}, has event starting at {} with an evdp of 0!\n'.format(self.station,self.BHN[0].stats.starttime))
            traveltime = model.get_travel_times(10,self.BHN[0].stats.sac.gcarc,[phase])[0].time
        else:
            tt = model.get_travel_times((self.BHN[0].stats.sac.evdp),self.BHN[0].stats.sac.gcarc,[phase])
            # print(self.BHN[0].stats.sac)
            # print(tt)
            try:
                traveltime = tt[0].time
            except IndexError:
                print('index Error')
                traveltime =None

        evt_time = obspy.UTCDateTime(year = self.BHN[0].stats.sac.nzyear, julday = self.BHN[0].stats.sac.nzjday,hour=self.BHN[0].stats.sac.nzhour,minute=self.BHN[0].stats.sac.nzmin,second=self.BHN[0].stats.sac.nzsec,microsecond=self.BHN[0].stats.sac.nzmsec)
        self.tt = evt_time + traveltime
        self.tt_rel = traveltime
        return traveltime

    def check_phase_dist(self,phase_to_check):
        """
        Function to test if the given phase is actually measureable!
        """
        if phase_to_check == 'SKS':
            if self.gcarc < 145.0:
                return True
            else:
                print('Event-Station distance further than 145 deg, too far for SKS')
                return False

        elif phase_to_check == 'SKKS':
            if self.gcarc >= 105.0:
                return True
            else:
                print('Event-Station distance less than 105 deg, too short for SKKS')
                return False
        elif phase_to_check == 'ScS':
            if self.gcarc <= 95.0:
                return True
            else:
                print('Event-Station distance greater than 105 deg, too far for ScS')
        else:
            print('Phase {} not SKS or SKKS'.format(phase_to_check))
            return False

    def process(self,synth=False,c1=0.01,c2=0.5,window=False):
        """
        Function to bandpass filter and trim the components
        Seismograms are trimmed so that they start 1 minute before the expected arrival and end 2 minutes after the arrival
        c1 - [Hz] Lower corner frequency
        c2 - [Hz] Upper corner frequency (N.B I have used c2 = 00.5 intially and am now trying c2 = 0.3 to see if that improves SNR without cutting oiut too mcuh singal)
        By default traces will be filtered between 0.01Hz-0.5Hz
        """
        # print(synth)
      # De-mean and detrend each component
        self.BHN.detrend(type='demean') #demeans the component
        self.BHE.detrend(type='demean')
        self.BHZ.detrend(type='demean')
        # print(self.BHN)
#       Detrend
        self.BHN.detrend(type='simple') #De-trends component
        self.BHE.detrend(type='simple') #De-trends component
        self.BHZ.detrend(type='simple') #De-trends component
#       Filter each component. Bandpass flag gives a bandpass-butterworth filter
        self.BHN.filter("bandpass",freqmin= c1, freqmax= c2,corners=2,zerophase=True)
        self.BHE.filter("bandpass",freqmin= c1, freqmax= c2,corners=2,zerophase=True)
        self.BHZ.filter("bandpass",freqmin= c1, freqmax= c2,corners=2,zerophase=True)
#       Now trim each component to the input length

        if synth == False: # We only need to trim and set the window length for real data, not synthetics made with sacsplitwave
#       Now set the trim

            print('Trim Traces')
            # print('tt = {}, Start = {}, End = {}'.format(self.tt,self.BHN[0].stats.starttime,self.BHN[0].stats.endtime ))
            t1 = (self.tt - 60) #I.e A minute before the arrival
            t2 = (self.tt+ 120) #I.e Two minutes after the arrival
            # print('t1 = {}, t2 = {}'.format(t1,t2))
            # self.BHN.trim(self.BHN[0].stats.starttime + t1,self.BHN[0].stats.starttime + t2)
            # self.BHE.trim(self.BHE[0].stats.starttime + t1,self.BHE[0].stats.starttime + t2)
            # self.BHZ.trim(self.BHZ[0].stats.starttime + t1,self.BHZ[0].stats.starttime + t2)
            self.BHN.trim(t1,t2)
            self.BHE.trim(t1,t2)
            self.BHZ.trim(t1,t2)
            # print('BX? records are assumed to already by short, so no trimming')
    #       Add windowing ranges to sac headers user0,user1,user2,user3 [start1,start2,end1,end2]
            # print(self.BHN[0].stats.sac.user0)
            if self.BHN[0].stats.sac.user0 == 0:
                print("Setting Window start/end ranges")
                # Set the range of window starttime (user0/user1)
                user0 = self.tt_rel - 15 #15 seconds before arrival
                user1 = self.tt_rel # t predicted arrival
        #       Set the raqnge of window endtime (user2/user3)
                user2 = self.tt_rel + 15 # 15 seconds after, gives a min window size of 20 seconds
                user3 = self.tt_rel + 30 # 30 seconds after, gives a max window size of 45 seconds
                keychain = {'user0':user0,'user1':user1,'user2':user2,'user3':user3}
                self.BHE[0].stats.sac.update(keychain)
                self.BHN[0].stats.sac.update(keychain)
                self.BHZ[0].stats.sac.update(keychain)
            else:
                print("Windows already set, user0-3 already set")
                print("User0 = ", self.BHN[0].stats.sac.user0)
                print("User1 = ",self.BHN[0].stats.sac.user1)
                print("User2 = ",self.BHN[0].stats.sac.user2)
                print("User3 = ",self.BHN[0].stats.sac.user3)

            if window == True:
                # Windowing code
                # Combine BHN and BHE to make a stream
                st = self.BHE + self.BHN
                # print(user0,user1,user2,user3)
                Windower = WindowPicker(st,user0,user1,user2,user3,self.tt_rel)
                # Windower.pick()
                if Windower.wbeg1 is None:
                    print("Skipping")
                    self.bad = True # Switch to tell sheba.py if we actually want to meausre this event
                else:
                    print("Windower Closed, adjusting window ranges")
                    (user0,user1,user2,user3) = Windower.wbeg1, Windower.wbeg2, Windower.wend1, Windower.wend2
                    self.bad = False
                # Set window ranges in SAC headers
                self.BHN[0].stats.sac.user0,self.BHN[0].stats.sac.user1,self.BHN[0].stats.sac.user2,self.BHN[0].stats.sac.user3 = (user0,user1,user2,user3)
                self.BHE[0].stats.sac.user0,self.BHE[0].stats.sac.user1,self.BHE[0].stats.sac.user2,self.BHE[0].stats.sac.user3 = (user0,user1,user2,user3)
                self.BHZ[0].stats.sac.user0,self.BHZ[0].stats.sac.user1,self.BHZ[0].stats.sac.user2,self.BHZ[0].stats.sac.user3 = (user0,user1,user2,user3)
        else:
            # print('Synthetics Used, Windows *should* be predefined')
            # print(self.BHN.stats.sac.user0,self.BHN.stats.sac.user1,self.BHN.stats.sac.user2,self.BHN.stats.sac.user3)
            pass

    def write_out(self,phase,label,path=None,synth=False):
        """
        Function to write the component seismograms to SAC files within the sheba directory structure so Sheba can access them
        station [str] - station code
        phase [str] - phase code for which seismic phase splitting is being meausured
        i [int] - counter for nukber of events (with the desired phase) at the station
        path [str] - path that you want the seismogrmas saved to.
        """
#       Now write out the three processed components
#       Naming depends on whether this is being executed as a test or within a loop
#       where a counter should be provided to prevent overwriting.
        if path is not None:
            self.BHN.write('{}/{}{}.BHN'.format(path,label,phase,self.ch),format='SAC',byteorder=1)
            self.BHE.write('{}/{}{}.BHE'.format(path,label,phase,self.ch),format='SAC',byteorder=1)
            self.BHZ.write('{}/{}{}.BHZ'.format(path,label,phase,self.ch),format='SAC',byteorder=1)
        elif synth == True:
            self.BHN.write('{}/{}SYNTH.BHN'.format(path,label,self.ch),format='SAC',byteorder=1)
            self.BHE.write('{}/{}SYNTH.BHE'.format(path,label,self.ch),format='SAC',byteorder=1)
            self.BHZ.write('{}/{}SYNTH.BHZ'.format(path,label,self.ch),format='SAC',byteorder=1)
        else:
            self.BHN.write('{}.BHN'.format(label,self.ch),format='SAC',byteorder=1)
            self.BHE.write('{}.BHE'.format(label,self.ch),format='SAC',byteorder=1)
            self.BHZ.write('{}.BHZ'.format(label,self.ch),format='SAC',byteorder=1)

    def plot_comp(self):
        """
        Quick Function to plot component together on one seismogram
        """
        st = self.BHN + self.BHE + self.BHZ
        st.plot(type='relative')

    # def gen_infile(self,path,label,phase,nwind=10,tlag_max=4.0):
    #
    #     os.chdir(path) # Make sure we are in the right directory
    #
    #     with open('sheba.in','w') as writer:
    #         writer.write('SHEBA.IN \n')
    #         writer.write('{}{} \n'.format(label,phase)) # write file prefix
    #         writer.write('{} \n'.format(self.BHE[0].stats.channel)) # Write channels for each component (E, N, Z order)
    #         writer.write('{} \n'.format(self.BHN[0].stats.channel))
    #         writer.write('{} \n'.format(self.BHZ[0].stats.channel))
    #         writer.write('1 \n') # Specifies Eigenvalue minimisation, replace with spol if transverse minimisation is desired (not supported here)
    #         writer.write('{:i} {:i} \n'.format(nwind,nwind))
    #         writer.write('{} \n'.format(tlag_max)) # sets max tlag in gridsearch
    #         writer.write('0')

    def sheba(self,phase,label,path=None,nwind=True,):
        """
        The big one! This function uses the subprocess module to host sac and then runs sheba as a SAC macro
        """
        print('Worker {} Passing {} into Sheba for {}. Time is {}'.format(current_process().pid,label,phase,time.ctime()))
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
            # print(s)
            out = p.communicate(s)
            # print(out[0])
        except CalledProcessError as err:
            print(err)