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
freq = [0.3,0.6]
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

finalevtlst = []
ccstackZ = []
Zstacklst = []
Rstacklst = []
trnum = []

for evtid in evtid_uni:
    evtid_rf = evtid[:-2]
    ievt = np.where(lst==evtid)[0][0]
    Cstalst = stalst[Cevtid==evtid]
    if len(Cstalst) < 4: continue
    ccstalst = []
    data = []
    for sta in Cstalst:
        if sta not in stalst_normal: continue
        try:
            tr = obspy.read(os.path.join(rotate_path,evtid_rf,'%s_%s.Z'%(evtid_rf,sta)))[0].copy()
        except:
            print("can't read %s, continue"%os.path.join(rotate_path,evtid,'%s_%s.Z'%(evtid,sta)))
            continue
        npts = tr.stats.npts
        ccstalst.append(sta)
        data.append(preprocess(tr,be=[-time_before,time_after],freq=wf[ievt,2:],win=wf[ievt,:2]))
    if len(ccstalst) < 3: continue
    sampling_rate = tr.stats.sampling_rate
    dist = tr.stats.sac.gcarc
    azi = tr.stats.sac.az
    srcdep = tr.stats.sac.evdp
    try: 
        F = radiationfactor(azi,srcdep,dist,fault[evtid[:-2]][0],fault[evtid[:-2]][1],fault[evtid[:-2]][2])
    except:
        print('no CMT for %s, continue'%evtid)
        continue
    time_axis = np.linspace(-time_before,time_after,num=npts,endpoint=True)
    data = np.array(data).T
    dcor,_,_,_,_ = cc.mccc(data,1./sampling_rate,ifwt=False)
    stackZ = np.zeros(npts)
    stackR = np.zeros(npts)
    for i in range(len(ccstalst)):
        sta = ccstalst[i]
        Z = obspy.read(os.path.join(rotate_path,evtid_rf,'%s_%s.Z'%(evtid_rf,sta)))[0].copy()
        R = obspy.read(os.path.join(rotate_path,evtid_rf,'%s_%s.R'%(evtid_rf,sta)))[0].copy()
        Z.filter('bandpass',freqmin=freq[0],freqmax=freq[1],zerophase=True)
        R.filter('bandpass',freqmin=freq[0],freqmax=freq[1],zerophase=True)
        stackZ += np.interp(time_axis,time_axis-dcor[i],Z.data)
        stackR += np.interp(time_axis,time_axis-dcor[i],R.data)
    stackZ /= (len(ccstalst)*F)
    stackR /= (len(ccstalst)*F)
    if isnormalize:
        Pnorm = np.abs(stackZ[np.where((time_axis>=ccwin[0])&(time_axis<=ccwin[1]))[0]]).max()
        stackZ = stackZ/Pnorm
        stackR = stackR/Pnorm
    finalevtlst.append(evtid)
    ccstackZ.append(stackZ[np.where((time_axis>=ccwin[0])&(time_axis<=ccwin[1]))[0]])
    Zstacklst.append(stackZ)
    Rstacklst.append(stackR)
    trnum.append(len(ccstalst))
    print(evtid+' succeeded')

time_axis = np.linspace(-time_before,time_after,num=npts,endpoint=True)
'''
stackZ = np.zeros(npts)
stackR = np.zeros(npts)
fig = plt.figure(figsize=(10,10))
ax1 = plt.axes([0.2,0.8,0.6,0.1])
ax2 = plt.axes([0.1,0.1,0.37,0.65])
ax3 = plt.axes([0.53,0.1,0.37,0.65])
ax1.set(xlim=[-1,10],xlabel='t[s]',title='overall PRF')
ax2.set(xlim=[-20,20],xlabel='t[s]',ylim=[-1,len(finalevtlst)+1],yticks=np.arange(0,len(finalevtlst)+1),title='Z')
ax3.set(xlim=[-20,20],xlabel='t[s]',ylim=[-1,len(finalevtlst)+1],yticks=np.arange(0,len(finalevtlst)+1),title='R')
'''
for align in ['none','maxamp']:#['none','maxamp','mccc']:
    stackZ = np.zeros(npts)
    stackR = np.zeros(npts)
    fig = plt.figure(figsize=(7,10))
    ax1 = plt.axes([0.2,0.8,0.6,0.1])
    ax2 = plt.axes([0.1,0.1,0.37,0.65])
    ax3 = plt.axes([0.53,0.1,0.37,0.65])
    ax1.set(xlim=[-1,10],xlabel='t[s]',title='overall PRF')
    ax2.set(xlim=[-20,20],xlabel='t[s]',ylim=[-1,len(finalevtlst)+1],yticks=np.arange(0,len(finalevtlst)+1),title='Z (%s)'%align)
    ax3.set(xlim=[-20,20],xlabel='t[s]',ylim=[-1,len(finalevtlst)+1],yticks=np.arange(0,len(finalevtlst)+1),title='R (%s)'%align)
    if align == 'mccc':
        dcor,_,_,_,_ = cc.mccc(np.array(ccstackZ).T,1./sampling_rate,ifwt=False)
        for i in range(len(finalevtlst)):
            ax2.plot(time_axis-dcor[i],Zstacklst[i]+i,'k-',linewidth=1)
            ax3.plot(time_axis-dcor[i],Rstacklst[i]+i,'k-',linewidth=1)
            stackZ += np.interp(time_axis,time_axis-dcor[i],Zstacklst[i])
            stackR += np.interp(time_axis,time_axis-dcor[i],Rstacklst[i])
    elif align == 'none':
        for i in range(len(finalevtlst)):
            ax2.plot(time_axis,Zstacklst[i]+i,'k-',linewidth=1)
            ax3.plot(time_axis,Rstacklst[i]+i,'k-',linewidth=1)
            stackZ += Zstacklst[i]
            stackR += Rstacklst[i]
    elif align == 'maxamp':
        for i in range(len(finalevtlst)):
            time_axis_part = time_axis[np.where((time_axis>=ccwin[0])&(time_axis<=ccwin[1]))[0]]
            stackZ_part = Zstacklst[i][np.where((time_axis>=ccwin[0])&(time_axis<=ccwin[1]))[0]]
            t_max = time_axis_part[stackZ_part.argmax()]
            ax2.plot(time_axis-t_max,Zstacklst[i]+i,'k-',linewidth=1)
            ax3.plot(time_axis-t_max,Rstacklst[i]+i,'k-',linewidth=1)
            stackZ += np.interp(time_axis,time_axis-t_max,Zstacklst[i])
            stackR += np.interp(time_axis,time_axis-t_max,Rstacklst[i])
    stackZ /= len(finalevtlst)
    stackR /= len(finalevtlst)
    ax2.plot(time_axis,stackZ+i+1,'r-',linewidth=2)
    ax3.plot(time_axis,stackR+i+1,'r-',linewidth=2)
    ax2.set(yticklabels=['']*len(finalevtlst)+['stack'])
    ax3.set(yticklabels=['']*len(finalevtlst)+['stack'])
    ax2.grid(); ax3.grid()
    rf_R = decov.decovt(stackZ,stackR,1./sampling_rate,time_beforeP,time_afterP,gauss)
    RF_time = np.linspace(-time_beforeP,time_afterP,num=len(rf_R),endpoint=True)
    ax1.plot(RF_time,rf_R,'k-',linewidth=2)
    ax1.grid()
    plt.savefig('test%.1f_%.1fhz_%s.png'%(freq[0],freq[1],align))
    plt.close('all')
