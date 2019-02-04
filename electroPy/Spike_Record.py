# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 10:04:17 2018

@author: ludovic.spaeth
"""
import neo 
import numpy as np
from scipy import signal
from matplotlib import pyplot as plt 

class Spike_Record:
    def __init__(self,filename,ext='winwcp',channel=0,leak_remove=True,plot=True):
        
        if ext == 'winwcp':
            r =neo.io.WinWcpIO(filename=filename)
            self.block = r.read_block(False,True) 
            
        else:
            print ('Other format than WinWcp are not implemented yet...')
            
    
    #def trace(self,channel,leak_remove=True,plot=True):  #Turns recording to a list sweep[[analog][time]]
        #fill with channel number
        #if leak_remove=False, trace will be loaded as it was recorded without leak compensation
        #kwarg for plot display or not 
        
        trace = []
            
        for i in range(len(self.block.segments)):
            
            matrix = np.zeros((2,len(np.ravel(self.block.segments[i].analogsignals[channel]))))
            
            if leak_remove==True:
                matrix[1,:]=np.squeeze(self.block.segments[i].analogsignals[channel])
                leak = np.mean(matrix[1,50:550])
                matrix[1,:] = matrix[1,:]-leak #trace without leak
            else:
                matrix[1,:]=np.squeeze(self.block.segments[i].analogsignals[channel]) #Trace with leak           
            
            matrix[0,:]=np.squeeze(self.block.segments[i].analogsignals[channel].times)#Time in seconds
        
            trace.append(matrix)
            
        self.trace=trace
        self.sampling_rate = self.block.segments[1].analogsignals[channel].sampling_rate
            
        if plot==True:
            fig,ax=plt.subplots(int(len(self.block.segments)/5),int(len(self.block.segments)/int(len(self.block.segments)/5)))
            plt.suptitle('Raw Recordings')
            fig.subplots_adjust(hspace=.5,wspace=.5)
            ax = ax.ravel()
            
            for i in range(len(self.block.segments)):
                ax[i].plot(trace[i][0],trace[i][1],linewidth=0.5)
                ax[i].set_title('Sweep #%s'%i)
                ax[i].set_ylabel('Amplitude')
                ax[i].set_xlabel('Time (s)')

            plt.show()
        
        #return trace
    
#---------------------------Now some more specific plots---------------------------------------------------------    
    
    def raster(self,threshold,distance=50,plot=True):
        #EXECUTE THE TRACE METHOD FIRST
        #Function to get spike time index  and raster plot
        #Threshold is the minimum spike size to be recognized
        #distance in points = minimum horizontal space between spikes 
       
        spike_times = []
        
        time = self.trace[0][0]  #report all sweeps on the first one for the time
        
        if plot == True:
            plt.figure()
            plt.suptitle('Raster Plot')
            plt.xlabel('Time (s)')
            plt.ylabel('Sweep #')
            plt.xlim(time[0],time[-1])
        
        for i in range(len(self.block.segments)):
            seg = self.trace[i][1] #the analogsignal
    
            spikes, _ = signal.find_peaks(seg,height=threshold,distance=distance) #Spike Indexes
            spike_idx = np.squeeze(spikes/self.sampling_rate) #Spikes times
            spike_times.append(spike_idx)
            
            if plot==True:
                plt.eventplot(spike_idx,orientation='horizontal',lineoffsets=-i,linelengths=0.8)
            
        return spike_times
    
    

            #plt.eventplot(spike_times,orientation='horizontal',lineoffsets=-i,linelengths=0.8)
            #plt.xlim(time[0],time[-1])
     
#------------------TEST CODE---------------------------------
'''       
filename = 'U:/RAW DATA/Controls/Granule Cell Spiking RT/GC_1.wcp'

neuron = Spike_Record(filename,'winwcp',0).trace(0)

from matplotlib import pyplot as plt

plt.plot(neuron[0][0],neuron[0][1])
'''








