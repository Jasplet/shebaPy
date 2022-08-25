import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.fftpack import fft
from scipy.signal import hilbert

# from .plots import plot_traces

class WindowPicker:
    """
    Picks a Window start/end range, for use with cluster analysis code
    """

    def __init__(self, st ,wbeg1 ,wbeg2 ,wend1 ,wend2 ,tt,**kwargs):
        '''
        Cosntructs ineractive window picker for a waveform
        Parameters
        ----------
        st : 
            obspy Stream object containing waveform data
        wbeg1 : float
            initial min value of window start range
        wbeg2 : float
            initial max value of window start range
        wend1 : float
            initial min value of window end range
        wend2 : float
            initial max value of window end range
        tt : float
            phase arrival time relative relative to event source time [s]

        Returns
        -------
        None.
        '''
        #t0 = 60 seconds before traveltime (aka the start of the trimming seismogram)
        self.st = st # Obspy stream containing BHN and BHE
        # st.plot()
        # if 't1' in st[0].stats.sac:
        #     self.tt = st[0].stats.sac['t1']
        # else:
        self.tt = tt
        t0 = tt - 60
        self.delta = st[0].stats.delta
        self.t = t0 + st[0].times()
        # make initial window ranges attributes
        (self.wbeg1,self.wbeg2,self.wend1,self.wend2) = (wbeg1,wbeg2,wend1,wend2)
        (self.x1,self.x2,self.x3,self.x4) = (wbeg1,wbeg2,wend1,wend2)
        # Base plot (before interactive stuff)
        self.fig = plt.figure(figsize = (10,8))
        yr = st[0].stats.sac['nzyear']
        jd = st[0].stats.sac['nzjday']
        hr = st[0].stats.sac['nzhour']
        mn = st[0].stats.sac['nzmin']
        sc = st[0].stats.sac['nzsec']
        ms = st[0].stats.sac['nzmsec']
        plt.suptitle(f'Station {st[0].stats.station}, Event Time {yr:4d}-{jd:03d} {hr:02d}:{mn:02d}:{sc:02d}.{ms:3d}')
        gs = gridspec.GridSpec(2,2)
        #self.ax1 = plt.subplot(gs[0,:]) # Top Row, for fft plot
        self.ax2 = plt.subplot(gs[0,:]) # Middle Row, for window picking
        self.ax3 = plt.subplot(gs[1,:]) # Bottom row, for envelopes
        #self.plot_fft()
        # Add seismograms
        self.ax2.plot(self.t, self.st[0].data,label=st[0].stats.channel, color='darkorange')
        self.ax2.plot(self.t, self.st[1].data,label=st[1].stats.channel, color='dodgerblue')
        self.ax3.set_xlabel('Time relative to origin (s)')
        # Add instantaneous amplitude envelopes to help pick out signal (should be envelope max at phase arrival)
        ht1 = hilbert(self.st[0].data)
        env1 = np.abs(ht1)
        ht2 = hilbert(self.st[1].data)
        env2 = np.abs(ht2)
        self.ax3.plot(self.t, env1, color='darkorange',linestyle='--', label=None)
        self.ax3.plot(self.t, env2, color='dodgerblue',linestyle='--', label=None)
        self.ax3.set_ylabel('Instantaneous amplitude')
        # Add legend
        self.ax2.legend()
        # window limit lines
        self.wbeg1line = self.ax2.axvline(self.wbeg1, linewidth=2, color='r', visible=True)
        self.wbeg2line = self.ax2.axvline(self.wbeg2, linewidth=2, color='r', visible=True)
        self.wend1line = self.ax2.axvline(self.wend1, linewidth=2, color='g', visible=True)
        self.wend2line = self.ax2.axvline(self.wend2, linewidth=2, color='g', visible=True)
        self.cursorline= self.ax2.axvline(100, linewidth=1, color='0.5', visible=False)

        self.pred_tt= self.ax2.axvline(self.tt, linewidth=1, color='k', visible=True)

        _, self.ydat = self.wbeg1line.get_data()
        
        # set limits
        self.lim_max = max([self.st[0].data.max(), self.st[1].data.max()]) * 1.1
        self.lim_min = min([self.st[0].data.min(), self.st[1].data.min()]) * 1.1
        # self.ax1.set_aspect('equal')
        self.ax2.set_ylim([self.lim_min,self.lim_max])
        self.ax2.set_xlim(t0,max(self.t)) # Set ylim in relative time (from stsrt of stream )
        self.ax3.set_xlim(t0,max(self.t))
        # Add some labels
        self.phaselabel = self.ax2.text(self.tt + 1,
                                        self.lim_max*0.8,"IASP91\nPred.\nArrival",
                                        multialignment='left')
        self.wbeg1label = self.ax2.text(self.wbeg1 - 3, self.lim_min*0.85, 'S', color='r', fontsize=14)
        self.wbeg2label = self.ax2.text(self.wbeg2 - 3, self.lim_min*0.85, 'F', color='r', fontsize=14)
        self.wend1label = self.ax2.text(self.wend1 - 3, self.lim_min*0.85, 'S', color='g', fontsize=14)
        self.wend2label = self.ax2.text(self.wend2 - 3, self.lim_min*0.85, 'F', color='g', fontsize=14)
        print("'a' & 'd' set the window beginnning range")
        print("'z' & 'c' set the window end range")
        self.connect()
        plt.tight_layout()
        plt.show()

    def plot_fft(self):
        '''
        Takes and fft of both components (within wbeg1 and wend2 and plots).
        This does not get updated by changes in windowing
        '''
        # Trim traces to intial wbeg and wend. Use delta to work out the indexes that correspond to the window start/end positions.
        st = self.st.copy()

        horz1 = st[0]
        horz2 = st[1]
        N = len(horz1) # number of samples in traces
        if N % 2 !=0:
            N = N -1

        print(N)
        # Set sample spacing in f-domain
        df = 1.0/(2.0*self.delta) # Sample frequency
        xf = np.linspace(0.0,df,int(N/2)) # Frequencies up to F_nyquist (N/2*df)
        # Take the fft
        fft_h1 = fft(horz1.data)
        fft_h2 = fft(horz2.data)
        #Now plot the spectra
        self.ax1.semilogy(xf[1:N//2], 2.0/N * np.abs(fft_h1[1:N//2]), color='darkorange')
        self.ax1.semilogy(xf[1:N//2], 2.0/N * np.abs(fft_h2[1:N//2]), color='dodgerblue')
        self.ax1.legend([st[0].stats.channel, st[1].stats.channel])
        self.ax1.set_xlabel('Frequency [Hz]')
        self.ax1.set_xlim([0.01, 0.5])
       
        # plt.show()

    def connect(self):
        '''
        Sets up interactive options in figure canvas
        '''
        # self.cidclick = self.fig.canvas.mpl_connect('button_press_event', self.click)
        self.cidmotion = self.fig.canvas.mpl_connect('motion_notify_event', self.motion)
        # self.cidrelease = self.canvas.mpl_connect('button_release_event', self.release)
        self.cidenter = self.fig.canvas.mpl_connect('axes_enter_event', self.enter)
        self.cidleave = self.fig.canvas.mpl_connect('axes_leave_event', self.leave)
        self.cidkey = self.fig.canvas.mpl_connect('key_press_event', self.keypress)

    def keypress(self,event):
        ''' Define a set of keypress responses
        'a' & 'd' set the window beginnning range
        'z' & 'c' set the window end range
        'q' exit the plot and returns the current WBEG, WEND
        The vertical line markers and annotations are redrawn after each Key Press
        
        '''
        if event.key == "a":
            print('WBEG Start')
            self.x1 = event.xdata
            print(self.ydat)
            self.wbeg1line.set_data(self.x1,self.ydat)
            self.wbeg1label.set_position((self.x1 - 3, self.lim_min*0.85))
            plt.draw()
            print(self.x1)
        elif event.key == "d":
            print('WBEG End')
            self.x2 = event.xdata
            self.wbeg2line.set_data(self.x2,self.ydat)
            self.wbeg2label.set_position((self.x2 - 3, self.lim_min*0.85))
            plt.draw()
            print(self.x2)
        elif event.key == "z":
            print('WEND Start')
            self.x3 = event.xdata
            self.wend1line.set_data(self.x3,self.ydat)
            self.wend1label.set_position((self.x3 - 3, self.lim_min*0.85))
            plt.draw()
            print(self.x3)
        elif event.key == "c":
            print('WEND End')
            self.x4 = event.xdata
            self.wend2line.set_data(self.x4,self.ydat)
            self.wend2label.set_position((self.x4 - 3, self.lim_min*0.85))
            plt.draw()
            print(self.x4)
        elif event.key == "w":
            print('Bad/noisey waveform, will tell SHEBA to skip')
            self.x1 = False
            self.disconnect()
        elif event.key == "q":
            print('DISCONNECT')
            self.disconnect()
        else:
            print(f'Key {event.key} unknown')

    def enter(self,event):
        if event.inaxes is not self.ax2: return
        x = event.xdata
        self.cursorline.set_data([x,x],self.ydat)
        self.cursorline.set_visible(True)
        plt.draw()

    def leave(self,event):
        if event.inaxes is not self.ax2: return
        self.cursorline.set_visible(False)
        plt.draw()

    def motion(self,event):
        if event.inaxes is not self.ax2: return
        x = event.xdata
        self.cursorline.set_data([x,x],self.ydat)
        plt.draw()

    def disconnect(self):
        '''
        disconnect all the stored connection ids
        '''
        # self.fig.canvas.mpl_disconnect(self.cidclick)
        self.fig.canvas.mpl_disconnect(self.cidmotion)
        self.fig.canvas.mpl_disconnect(self.cidenter)
        self.fig.canvas.mpl_disconnect(self.cidleave)
        if self.x1 == False:
            print('Bad waveform that we want to skip')
            self.wbeg1, self.wbeg2, self.wend1,self.wend2 = (None,None,None,None)
        else:
            self.wbeg1, self.wbeg2, self.wend1,self.wend2 = sorted((self.x1, self.x2,self.x3,self.x4))
        plt.close()
