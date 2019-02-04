# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 17:52:37 2018

@author: Sam Garcia (adapted by ludovic.spaeth)
"""

from numpy import inf
import numpy 
from scipy import *
from electroPy import AnalogSignal
from scipy.signal import resample
from scipy.fftpack import fft,ifft,fftshift

# global for caching wf
global cache_for_wf
cache_for_wf = None
global signature_for_wf
signature_for_wf = ''

if True:

    def generate_wavelet_fourier(len_wavelet,
                                 f_start,
                                 f_stop,
                                 deltafreq,
                                 sampling_rate,
                                 f0,
                                 normalisation):
        
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
                wf = self.generate_wavelet_fourier(sig.size,max(f_start,deltafreq),min(f_stop,ana.sampling_rate/2.),deltafreq,sampling_rate,f0,normalisation)
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



class TimeFreq():
    doc2="""
    *TimeFreq*
    
    """
    docparam = """
    
    Params:
     :f_start: lkjlkj
    
    """
    
    __doc__ = doc2+docparam
    
    def __init__(self,
                        anaSig,
                        method = 'convolution_freq',
                        f_start=5.,
                        f_stop=20.,
                        deltafreq = 0.5,
                        sampling_rate = 10000.,
                        t_start = -inf, 
                        t_stop = inf,
                        f0=2.5, 
                        normalisation = 0.,
                        **kargs
                        ):
                        
        self.anaSig = anaSig
        self.method = method
        self.f_start=f_start
        self.f_stop=f_stop
        self.deltafreq = deltafreq
        self.sampling_rate = sampling_rate
        self.t_start = t_start
        self.t_stop = t_stop
        self.f0=f0
        self.normalisation = normalisation
        
        self.t_start = max(self.t_start , self.anaSig.t_start)
        self.t_stop = min(self.t_stop , self.anaSig.t()[-1]+1./self.anaSig.sampling_rate )
        
        self._map = None
        self._t = None
        self._f = None
        
        if self.method == 'convolution_freq':
            self._wf = None
            self.subAnaSig = None
    def compute_time_vector(self) :
        return numpy.arange(len(self.subAnaSig.signal), dtype = 'f8')/self.sampling_rate + self.t_start
    def compute_freq_vector(self) :
        return numpy.arange(self.f_start,self.f_stop,self.deltafreq, dtype = 'f8')
    def t(self):
        if self._t==None:
            self._t=self.compute_time_vector()
        return self._t
    
    def f(self):
        if self._f==None:
            self._f=self.compute_freq_vector()
        return self._f
    
        

    
    
    
    @property
    def map(self):
        if self._map is None:
            self.recomputeMap()
        return self._map
    
    
    
    
    def recomputeMap(self):
        """
        Compute or recompute a map
        """
        if self.subAnaSig is None:
        #~ if True:
            sig=self.anaSig.signal[(self.anaSig.t()>=self.t_start)&(self.anaSig.t()<self.t_stop)]
            if self.sampling_rate != self.anaSig.sampling_rate :
                sig=resample(sig,sig.size*self.sampling_rate/self.anaSig.sampling_rate)
            self.subAnaSig = AnalogSignal( signal = sig,
                                                            sampling_rate = self.sampling_rate,
                                                            t_start = self.t_start,
                                                        )
        
        if self.method == 'convolution_freq':
            if self._wf is None :
                global signature_for_wf
                global cache_for_wf
                signature = '%d %f %f %f %f %f %f' % (self.subAnaSig.signal.size,
                                                                                self.f_start,
                                                                                self.f_stop,
                                                                                self.deltafreq,
                                                                                self.sampling_rate,
                                                                                self.f0,
                                                                                self.normalisation,
                                                                                )
                if signature != signature_for_wf:
                    cache_for_wf= generate_wavelet_fourier(len_wavelet=self.subAnaSig.signal.size,
                                                                                f_start=self.f_start,
                                                                                f_stop=self.f_stop,
                                                                                deltafreq=self.deltafreq,
                                                                                sampling_rate=self.sampling_rate,
                                                                                f0=self.f0,
                                                                                normalisation = self.normalisation
                                                                                )
                    signature_for_wf = signature
                
                self._wf = cache_for_wf
            
            
            self._map = compute_morlet_scalogram(self.subAnaSig,wf = self._wf )
            
        
    def plotMap(self, ax,
                                    colorbar = True,
                                    cax =None,
                                    orientation='horizontal',
                                    **kargs):
        """
        
        ax : a matplotlib axes
        
        """
        im = ax.imshow(abs(self.map).transpose(),
                                    interpolation='nearest', 
                                    extent=(self.t_start, self.t_stop, self.f_start-self.deltafreq/2., self.f_stop-self.deltafreq/2.),
                                    #origin ='lower' ,
                                    )
        if colorbar:
            if cax is None:
                ax.figure.colorbar(im)
            else:
                ax.figure.colorbar(im,ax = ax, cax = cax ,orientation=orientation)
            
                
        return im
    