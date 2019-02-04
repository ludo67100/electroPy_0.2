# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 16:17:24 2018

@author: Sam Garcia (adapted by ludovic.spaeth)

Script for convolution in scalogram design 

"""

from scipy import *
from scipy.signal import resample
from scipy.fftpack import fft,ifft,fftshift


class Convolution():
    
    def __init__(self):
        print ('Convolution module loaded')
        

        
    def generate_wavelet_fourier(len_wavelet,
                                 f_start,
                                 f_stop,
                                 deltafreq,
                                 sampling_rate,
                                 f0,
                                 normalization):
        
        '''
        Computes the wavelet coefficients at all scales and makes corresponding Fourier transform
        When different signal scalograms are computed with the exact same coefficient
        this function can be excecuted once and the result is passed directly to compute_morlet_scalogram
        
        Output = wf : Fourier transform of the waverlet coefficients (after weighting), Fourier freq at first
        '''
        
        #Computes final map scales
        scales = f0/arange(f_start,f_stop,deltafreq)*sampling_rate
       
        #compute wavelet coeffs at all scales
        xi = arange(-len_wavelet/2.0,len_wavelet/2.0)
        xsd = xi[:,newaxis] / scales
        wavelet_coefs = exp(complex(1j)*2.0*pi*f0*xsd)*exp(-power(xsd,2)/2.)
        
        weighting_function = lambda x: x**(-(1.0))
        wavelet_coefs = wavelet_coefs*weighting_function(scales[newaxis,:])
        
        #Transform the wavelet into the Fourier domain 
        wf=fft(wavelet_coefs,axis=0)
        wf=wf.conj() #Used to be an error here
        
        
        return wf 
    
    
    def compute_morlet_scalogram(ana,
                                 f_start=5,
                                 f_stop=100,
                                 deltafreq=1,
                                 sampling_rate=200,
                                 t_start=-inf,
                                 t_stop=inf,
                                 f0=2.5,
                                 normalisation=0,
                                 wf=None):
        
        #Reduce the signal to the limits
        sig = ana.signal[(ana.t()>=t_start)&(ana.t()<=t_stop)]
        
        if sig.size>0:
            if wf is None:
                if ana.sampling_rate != sampling_rate:
                    sig = resample(sig,sig.size*sampling_rate/ana.sampling_rate)
                wf = Convolution.generate_wavelet_fourier(sig.size,max(f_start,deltafreq),min(f_stop,ana.sampling_rate/2.),deltafreq,sampling_rate,f0,normalisation)
            else:
                if sig.size != wf.shape[0]:
                    sig=resample(sig,wf.shape[0])
                    
            #Transform the signal in Fourier domain 
            sigf=fft(sig)
            
            #Convolve (mult. in Fourier space)
            wt_tmp = ifft(sigf[:,newaxis]*wf,axis=0)
            
            #Shift output from ifft
            wt = fftshift(wt_tmp,axes=[0])
            
        else: 
            scales = f0/arange(f_start,f_stop,deltafreq)*sampling_rate
            wt = empty((0,scales.size),dtype='complex')
            
        return wt 
                    