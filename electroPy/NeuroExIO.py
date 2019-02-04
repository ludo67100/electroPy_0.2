# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 17:14:53 2018

Class to load .nex files as np matrices. We consider raw signals (ch 0 to 15), spike filtered (16 to 31) and
LFP signal (ch 32 to 48)

returs info and data such as sampling rate, time vector, nb of channels and data itself 

@author: ludovic.spaeth
"""


class NeuroExIO:
    
    def __init__(self,filename):
        from neo.io import NeuroExplorerIO as NeuX
        r = NeuX(filename)
        self.block = r.read_segment()
        self.channels = len(self.block.analogsignals) #How many channels
        self.raw_points = len(self.block.analogsignals[0]) #how many points in the raw traces
        self.filt_points = len(self.block.analogsignals[17]) #how many points in the raw traces
        self.lfp_points = len(self.block.analogsignals[33]) #how many points in the raw traces
        
        
        
    def sampling_rate(self): #Returns sampling rate
        
        raw_sr = float(self.block.analogsignals[0].sampling_rate)
        filt_sr = float(self.block.analogsignals[17].sampling_rate)
        lfp_sr = float(self.block.analogsignals[33].sampling_rate)
        
        return [raw_sr,filt_sr,lfp_sr]
        
    
    def channels(self): #Returns total amout of channels in file, should be 3*n_probe_channels
        return int(self.channels)
    
    
    def time(self): #Returns time segment for data           
        import numpy as np        
        raw_time = np.ravel(self.block.analogsignals[0].times) #All ch should have same time  
        filt_time = np.ravel(self.block.analogsignals[17].times) #All ch should have same time  
        lfp_time = np.ravel(self.block.analogsignals[33].times) #All ch should have same time  
        return [raw_time,filt_time,lfp_time]
    
    
    def units(self): #Returns units [time,data] as str
        time_unit = self.block.analogsignals[0].times.units
        trace_unit = self.block.analogsignals[0].units
        return [str(time_unit)[-1:],str(trace_unit)[-2:]]        
    
    
    def raw_record(self): #Returns raw records as matrix (channel,data)
        import numpy as np
        raw_records = np.zeros((int(self.channels/3.0),self.raw_points))
        
        raw_indexes = np.arange(0,16,1) #Raw records are stored on ch 0 to 15 in the global file
        
        for i in raw_indexes:
            raw_records[i,:] = np.ravel(np.squeeze(self.block.analogsignals[raw_indexes[i]]))
            
        return raw_records
    
    def filt_record(self): #Returns spike-filtered data
        import numpy as np 
        filt_records = np.zeros((int(self.channels/3.0),self.filt_points))
        
        filt_indexes = np.arange(16,16+16,1) #Channels 16 to 31

        for i in range(16):
            filt_records[i,:] = np.ravel(np.squeeze(self.block.analogsignals[filt_indexes[i]]))
            
        return filt_records
       
        
    def LFP_record(self): #Returns LFP filtered data
        import numpy as np
        LFP_records = np.zeros((int(self.channels/3.0),self.lfp_points))
        
        LFP_indexes = np.arange(32,32+16,1) #Channels  32 to 48
        
        for i in range(16):
            LFP_records[i,:] = np.ravel(np.squeeze(self.block.analogsignals[LFP_indexes[i]]))
            
        return LFP_records
            
            
            
            
            
#----------TRIAL----------------------------------------------
 
filepath = r"H:\Federica\2018-08-10T13-54-41McsRecording.nex"

data = NeuroExIO(filepath)

time = data.time()

raw = data.raw_record()

filt = data.filt_record()

LFP = data.LFP_record()

from matplotlib import pyplot as plt

fig,ax = plt.subplots(6,3,sharex=True,sharey=True)

i = 0
for ax in ax.reshape(-1):

    
    
    ax.plot(time[0],raw[i],'r-')
    
    i = i+1
    
    
