# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 16:48:23 2018

This module is used to compute spike detetion with electroPy GUI
Input should be self.tagged_traces from GUI window, I should add an assertion error to check this 


@author: ludovic.spaeth
"""

'''
#fake data & running example

import neo
import scipy as sp
import quantities as pq

#fake analogsignal
AnalogSignal = neo.core.AnalogSignal(sp.randn(10000)*pq.uV,t_start=0*pq.s,sampling_rate=25*pq.kHz)

#Matrix same shape as what gets out of the Gui (nsweeps,npoints,2)
import numpy as np 
data = np.zeros((5,10000,2))

for i in range(5):
    data[i,:,1]=np.ravel(neo.core.AnalogSignal(sp.randn(10000)*pq.uV,t_start=0*pq.s,sampling_rate=25*pq.kHz))

from matplotlib import pyplot as plt 

plt.plot(data[1,:,1])
plt.show()
'''

class SpikeDetection:
    
    def __init__(self,matrix,leak_remove=True):
        print ('Spike Detect Module Loaded')
           
        self.data = matrix #Container for the raw data
        
    def spiketimes(self,threshold,distance=10):   
            
        import numpy as np
        import scipy.signal as sp 
        
        self.spike_times,self.sweep,self.threshold = [],[],[] #to return everything
        
        for i in range(self.data.shape[0]):
            seg = self.data[i,:,1] #the analogsignal
    
            spikes, _ = sp.find_peaks(seg,height=threshold,distance=distance) #Spike Indexes
            spike_idx = np.squeeze(spikes) #Spikes times
            self.spike_times.append(spike_idx)
            self.sweep.append(np.ones(len(spike_idx))*i)
            self.threshold.append(np.ones(len(spike_idx))*threshold)
            
        return (np.ravel(np.asarray(self.spike_times)),np.ravel(np.asarray(self.threshold)),np.ravel(np.asarray(self.sweep)))
            
'''
#Running example
a = SpikeDetection(data).spiketimes(threshold=2.5)   

plt.scatter(a[0][1],a[1][1],color='orange')     
'''


        