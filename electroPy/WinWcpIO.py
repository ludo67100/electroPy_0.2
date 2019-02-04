# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 11:51:49 2018

This class loads Winwcp files for GUI usage
Allows to load full file (init) as block (with quantities inside) 
Or a single sweep at a time in ndarray format 
Or full record as ndarray format, cause, let's be honnest, it's better to work with arrays in Python 

@author: ludovic.spaeth
"""

#filename = 'U:/RAW DATA/Mapping PHD/AD-02-05-2016/160502_001.wcp'

class WinWcpIO:
    
    def __init__(self,filename):
        from neo.io import WinWcpIO as win
        r = win(filename)
        self.bl = r.read_block(False,True)
        self.n = len(self.bl.segments) #number of sweeps
        self.points = len(self.bl.segments[0].analogsignals[0]) #Points in one sweep
    
    #Method to get sampling rate    
    def sampling_rate(self,sweep,channel):
        return float(self.bl.segments[sweep].analogsignals[channel].sampling_rate)
                
    #Method to get one sweep at a time, time in seconds [0] and trace [1] 
    def sweep(self,sweep,channel):   
        import numpy as np 
        trace = self.bl.segments[sweep].analogsignals[channel]
        time = trace.times        
        return np.ravel(np.array(time)),np.ravel(np.array(trace))
    
    #Method to get whole file as a np.matrix
    def whole_file(self,channel):
        import numpy as np 
        print ('Whole file has beel loaded succesfully')
        _analogsignals = np.zeros((self.n,self.points,2)) #z0 for time and z1 for current
        
        for i in range(self.n):
            _analogsignals[i,:,0] = np.ravel(self.bl.segments[i].analogsignals[channel].times) #Timepoints
            _analogsignals[i,:,1] = np.ravel(self.bl.segments[i].analogsignals[channel]) #Current points
            
        return _analogsignals
    

    
