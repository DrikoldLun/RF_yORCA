#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 19:40:27 2021

@author: lun
"""

import numpy as np
from radiation import radiationfactor
import configparser, os, sys
import obspy
sys.path.append('../RFscript')
import decov, geo
import cc
import matplotlib.pylab as plt

def preprocess(tr,be=[-40,40],freq=[0.1,2],win=[-3,5]):
    npts = tr.stats.npts
    time = np.linspace(be[0],be[1],num=npts,endpoint=True)
    tr.filter('bandpass',freqmin=freq[0],freqmax=freq[1],zerophase=True)
    tr.data = tr.data[np.where((time>=win[0])&(time<=win[1]))]
    tr.taper(max_percentage=0.05,type='hann')
    return tr.data

stalst_noise = ['CC06','EE04','CC11','WC05','WW03'] #0.8hz noise
stalst_normal = ['CC02','CC04','CC05','CC07','CC08','WC03','WC04','WC02'] #normal stations
#stalst_normal += stalst_noise
mag = [6,10]
ccwin = [-3,3]
freq = [0.1,0.5]
parafile = '../paraRF.cfg'
CMTdat = 'CMT_%.1f_%.1f.dat'%(mag[0],mag[1])
evtlog_file = 'evt_Z.log'
dat = 'Z_mag_%.1f_%.1f.dat'%(mag[0],mag[1])
wf_file = 'Z_mag_%.1f_%.1f_winfreq.dat'%(mag[0],mag[1])
isnormalize = True
#align = 'mccc'

Cevtid = np.loadtxt(dat,usecols=0,dtype='str')
wf = np.loadtxt(wf_file)
evtid_uni = np.unique(Cevtid)
stalst = np.loadtxt(dat,usecols=4,dtype='str')
lst = np.loadtxt(evtlog_file,usecols=0,dtype='str')
maglog = np.loadtxt(evtlog_file,usecols=4)
lst = lst[(maglog>=mag[0])&(maglog<=mag[1])]
fault = {}
with open(CMTdat,'r') as f:
    for line in f.readlines():
        line = line.replace('\n','').strip(' ').split()
        fault[line[-1][:-1]] = [float(tmp) for tmp in line[2:5]]
    f.close()

############################
# ------- Get Para ------- #
############################
config = configparser.ConfigParser()
prefix = parafile.replace(os.path.basename(parafile),'')
if os.path.isfile(parafile):
    config.read(parafile)
else:
    print("head file doesn't exist!")
out_path = os.path.join(prefix,config.get('path','out_path'))
rotate_path = config.get('path','rotate_dir')+'_event'
rotate_path = os.path.join(out_path,rotate_path)
time_beforeP = config.getfloat('para','time_beforeP')
time_afterP = config.getfloat('para','time_afterP')
time_before = config.getfloat('para','RF_calcwinlen')
time_after = config.getfloat('para','RF_calcwinlen')
gauss = config.getfloat('para','gauss')
#freqmin = config.getfloat('para','freqmin')
#freqmax = config.getfloat('para','freqmax')
#dismin = config.getfloat('para','dismin')
#dismax = config.getfloat('para','dismax')
#magmin = config.getfloat('para','magmin')
#samprate = config.getfloat('para','samprate')
#noisegate = config.getfloat('para','noisegate')
#RF_calcwinl = config.getfloat('para','RF_calcwinl')
#RF_calcwinr = config.getfloat('para','RF_calcwinr')

Zlst = []
Rlst = []

for evtid in evtid_uni:
    evtid_rf = evtid[:-2]
    Cstalst = stalst[Cevtid==evtid]
    try:
        tr = obspy.read(os.path.join(rotate_path,evtid_rf,'%s_*.Z'%(evtid_rf)))[0]
    except:
        print("can't read %s, continue"%evtid)
        continue
    sampling_rate = tr.stats.sampling_rate
    dist = tr.stats.sac.gcarc
    azi = tr.stats.sac.az
    srcdep = tr.stats.sac.evdp
    npts = tr.stats.npts
    time_axis = np.linspace(-time_before,time_after,num=npts,endpoint=True)
    try:
        F = radiationfactor(azi,srcdep,dist,fault[evtid[:-2]][0],fault[evtid[:-2]][1],fault[evtid[:-2]][2])
    except:
        print('no CMT for %s, continue'%evtid)
        continue    

    for sta in Cstalst:
        if sta not in stalst_normal: continue
        try:
            Z = obspy.read(os.path.join(rotate_path,evtid_rf,'%s_%s.Z'%(evtid_rf,sta)))[0].copy()
            R = obspy.read(os.path.join(rotate_path,evtid_rf,'%s_%s.R'%(evtid_rf,sta)))[0].copy()
        except:
            print("can't read %s, continue"%os.path.join(rotate_path,evtid,'%s_%s.Z'%(evtid,sta)))
            continue
        Z.filter('bandpass',freqmin=freq[0],freqmax=freq[1],zerophase=True)
        R.filter('bandpass',freqmin=freq[0],freqmax=freq[1],zerophase=True)
        Z.data /= F
        R.data /= F
        if isnormalize:
            Pnorm = np.abs(Z.data[np.where((time_axis>=ccwin[0])&(time_axis<=ccwin[1]))[0]]).max()
            Z.data = Z.data/Pnorm
            R.data = R.data/Pnorm
        Zlst.append(Z.data)
        Rlst.append(R.data)
    print(evtid+' succeeded')

for align in ['none','maxamp']:
    stackZ = np.zeros(npts)
    stackR = np.zeros(npts)
    fig = plt.figure(figsize=(7,6))
    ax1 = plt.axes([0.1,0.55,0.8,0.35])
    ax2 = plt.axes([0.1,0.1,0.8,0.35])
    ax1.set(xlim=[-1,10],xlabel='t[s]')
    ax2.set(xlim=[-20,20],xlabel='t[s]')
    if align == 'none':
        for i in range(len(Zlst)):
            stackZ += Zlst[i]
            stackR += Rlst[i]
    elif align == 'maxamp':
        for i in range(len(Zlst)):
            time_axis_part = time_axis[np.where((time_axis>=ccwin[0])&(time_axis<=ccwin[1]))[0]]
            Z_part = Zlst[i][np.where((time_axis>=ccwin[0])&(time_axis<=ccwin[1]))[0]]
            t_max = time_axis_part[Z_part.argmax()]
            stackZ += np.interp(time_axis,time_axis-t_max,Zlst[i])
            stackR += np.interp(time_axis,time_axis-t_max,Rlst[i])
    stackZ /= len(Zlst)
    stackR /= len(Zlst)
    rf_R = decov.decovt(stackZ,stackR,1./sampling_rate,time_beforeP,time_afterP,gauss)
    RF_time = np.linspace(-time_beforeP,time_afterP,num=len(rf_R),endpoint=True)
    ax1.plot(RF_time,rf_R,'k-',linewidth=2,label='RF')
    ax2.plot(time_axis,stackZ,'r-',linewidth=2,label='Z')
    ax2.plot(time_axis,stackR,'b-',linewidth=2,label='R')
    ax1.grid(); ax1.legend(loc='upper right')
    ax2.grid(); ax2.legend(loc='upper right')
    ax2.set_xlim([-1,10])
    plt.suptitle('%.1fto%.1f[Hz] trnum:%d align:%s'%(freq[0],freq[1],len(Zlst),align))
    plt.savefig('directstack%.1f_%.1fhz_%s.png'%(freq[0],freq[1],align))
    plt.close('all')
