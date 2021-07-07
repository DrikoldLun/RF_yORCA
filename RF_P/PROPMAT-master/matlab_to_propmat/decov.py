#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 10:53:49 2019

@author: lun
"""

import numpy as np
from scipy.io import loadmat
from scipy.fftpack import fft, ifft, ifftshift
from obspy.signal.util import next_pow_2
import matplotlib.pylab as plt

def gauss_spectrum(dt,nt,alpha):
    #alpha: gaussian para
    df = 1./(dt*nt)
    nf_mid = int(np.ceil(0.5*nt))
    w = 2*np.pi*df*np.arange(0,nf_mid)
    
    gauss = np.zeros(nt,dtype='float')
    gauss[:nf_mid] = np.exp(-np.power(0.5*w/alpha,2))/dt
    gauss[-nf_mid+1:] = np.flipud(gauss[1:nf_mid])
    #gauss /= np.max(ifft(gauss))
    
    return gauss

def decovt(P,D,dt,time_beforeP=10,time_afterP=30,alpha=2,accept_mis=1e-8,itmax=1e3,isplot=False):
    #alpha: gaussian para
    P = np.array(P)
    D = np.array(D)
    if P.ndim != 1 or D.ndim != 1:
        raise ValueError('Only 1 dimension arrays are supported.')
    if len(P) != len(D):
        raise ValueError('The length of Parent and Daughter is not equal.')
        
    misfit_old = 9999999999999999;
    misfit = np.sqrt(np.power(P-D,2).sum())
    
    RF = np.zeros(len(P)*2-1)
    D_cur = D.copy()
    itnum = 0
    
    #iteration deconvolution
    auto_corr = np.correlate(P,P,'full')
    while itnum <= itmax and misfit_old > accept_mis*auto_corr.max():
        amp_corr = np.correlate(D_cur,P,'full')
        ind = np.abs(amp_corr).argmax()
        amp_rf = amp_corr[ind]/auto_corr[len(auto_corr)//2]
        RF[ind] += amp_rf
        D_sub = np.convolve(P,RF,'valid')
        D_cur = D - D_sub
        misfit_old = misfit
        misfit = np.sqrt(np.power(D_cur,2).sum())/len(D_cur)
        itnum += 1
    
    RF_Time = (np.arange(0,len(amp_corr))-(len(D_cur)+len(P)-1)//2)*dt
    RF[np.where((RF_Time<-time_beforeP)|(RF_Time>time_afterP))[0]] = 0
    
   
    '''
    #shift RF (set max amplitude as the zero point)
    O = np.argwhere(RF_Time==0)[0][0]
    nb, ne = RF.argmax()-O, len(RF)+RF.argmax()-O
    if nb < 0: nb = 0
    else: ne = len(RF)
    modified_RF = np.zeros(len(RF))
    modified_RF[0:ne-nb] = RF[nb:ne]
    RF = modified_RF
    '''
    
    #cut RF
    RF = RF[np.where((RF_Time>=-time_beforeP)&(RF_Time<=time_afterP))[0]]
    RF_Time = RF_Time[np.where((RF_Time>=-time_beforeP)&(RF_Time<=time_afterP))[0]]
    
    #gauss pulse convolution
    gauss = gauss_spectrum(dt,len(RF),alpha)
    RF = ifft(fft(RF)*gauss).real
    
    return RF

def decovf(P,D,dt,t_P,time_beforeP=10,time_afterP=30,alpha=2,c=0.1,isplot=False):
    P = np.array(P)
    D = np.array(D)
    if P.ndim != 1 or D.ndim != 1:
        raise ValueError('Only 1 dimension arrays are supported.')
    if len(P) != len(D):
        raise ValueError('The length of Parent and Daughter is not equal.')
    
    # initiation
    nfft = next_pow_2(len(P))
    P_tmp = np.zeros(nfft)
    D_tmp = np.zeros(nfft)
    P_tmp[:len(P)] = P
    D_tmp[:len(P)] = D
    gauss = gauss_spectrum(dt, nfft, alpha)
    
    # calculate frequency spectrum
    Pf = fft(P_tmp)
    Df = fft(D_tmp)
    
    #calculate denominator
    Pf_2 = Pf*Pf.conjugate()
    
    # add water level
    wl = np.max(Pf_2.real)*c
    Pf_2[np.where(Pf_2<wl)[0]] = wl
    lower = Pf_2
    
    # numerator
    upper = Df*Pf.conjugate()*gauss
    
    #calculate RF
    RFf = upper/lower
    RF = ifftshift(ifft(RFf,nfft).real)
    if t_P > 0:
        Pn_left, Pn_right = int((t_P-5)/dt-1), int((t_P+15)/dt)
        O = np.arange(Pn_left,Pn_right)[RF[Pn_left:Pn_right].argmax()]
    else:
        O = RF.argmax()
    RF_b, RF_e = int(O-time_beforeP/dt), int(O+time_afterP/dt)
    RF = RF[RF_b:RF_e+1]
    
    return RF

def test():
    path = './RFdata.mat'
    struct = loadmat(path)['data']['datZRT'][0][0]
    P = np.array(struct[:,0]).reshape(1,-1)[0]
    D = np.array(struct[:,1]).reshape(1,-1)[0]
    RF_t = decovt(P,D,0.1,10,30,2)
    RF_f = decovf(P,D,0.1,0,10,30,2)
    RF_time = np.arange(-10,30.1,0.1)
    fig = plt.figure(figsize=(10,3))
    plt.plot(RF_time,RF_t,'r-',label='decovt')
    plt.plot(RF_time,RF_f,'b-',label='decovf')
    plt.legend(loc='upper right')
    plt.xlabel('t[s]')
    plt.ylabel('Amp')
    plt.grid()
    plt.show()
