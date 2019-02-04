# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 15:28:46 2018

@author: Sam Garcia (adapted by ludovic.spaeth)

Modified AnalogSignal class from OpenElectrophy to create AnalogSignal object from data 

"""

import numpy

class AnalogSignal(object):
    
    def __init__(self,signal=None,time_vector=None,channel=None,name=None,
                 sampling_rate=1.0,t_start=0,t_stop=None,dt=None,**kargs):
        
        self.signal = signal
        self.channel = channel
        self.name = name
        self.sampling_rate = float(sampling_rate)
        self.t_start = float(t_start)
        self.t_stop = t_stop
        
        #Default for signal is empty array
        if self.signal is None:
            self.signal = numpy.array([])
            
        #Override sampling rate if dt is specified
        if dt is not None:
            self.sampling_rate = 1.0 / dt
        if self.sampling_rate == 0:
            raise(ValueError('Sampling rate cannot be equal to zero'))
            
        #Calculate self.t_stop
        if self.t_stop is None:
            self.t_stop = self.t_start + len(self.signal)/self.sampling_rate
            
        #Initialize variable for time array (change from original script)
        self._t = time_vector
        
    def compute_time_vector(self):
        return numpy.arange(len(self.signal),dtype='f8')/self.sampling_rate + self.t_start
    
    def t(self):
        self._t=self.compute_time_vector()
        return self._t

    def max(self):
        return self.signal.max()
    
    def min(self):
        return self.signal.min()
    
    def mean(self):
        return numpy.mean(self.signal)
    
    def time_slice(self,t_start,t_stop):
        #Get time axis and also trigger recompute if needed
        t = self.t()
        
        assert t_stop > t_start, "t_stop must be > t_start"
        
        #Find the int indices of self.signal closest to resquested limits
        i_start = int(numpy.rint((t_start - self.t_start) * self.sampling_rate))
        
        if i_start < 0 :
            print ("Warning: you requested data before signal starts")
            
        #Add one so that it is inclusive of t_stop
        i_stop = int(numpy.rint((t_stop - self.t_start) * self.sampling_rate)) + 1
        if i_stop > len(self.signal):
            print ("Warning : you requested data after the signal ended")
            i_stop = len(self.signal)
            
        #Slice the signal
        signal = self.signal[i_start:i_stop]
        
        #Create a new AnalogSignal with the specified data and the correct underlying time axis 
        result = AnalogSignal(signal=signal,sampling_rate=self.sampling_rate,t_start=t[i_start])
        
        return result 
            
            
#'''--------------------TEST----(should be commented)-----------------------'''
#
#from electroPy import HdF5IO
#
##Where is the file ? The file is in the kitchen
#path = 'H:/Federica/2018-08-10T13-54-41McsRecording.h5'
#
##Opening .h5 file
#data = HdF5IO(path) 
#
##Define channel of interest, t_start and time_vector and also the sampling rate
#channel_0 = data.raw_record()[0]
#t_start = data.raw_time()[0]
#time_vector = data.raw_time
#sampling_rate = data.raw_sampling_rate()
#
##Create analogsignal from data
#anasig = AnalogSignal(signal=channel_0,sampling_rate=sampling_rate,t_start=t_start,time_vector=time_vector)
#
##Show me this fucking analog signal 
#from matplotlib import pyplot as plt 
#plt.figure()
#plt.plot(anasig.t(),anasig.signal)
        
        
        