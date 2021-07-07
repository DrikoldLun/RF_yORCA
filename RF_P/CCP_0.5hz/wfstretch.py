#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 26 09:08:44 2020

@author: lun
"""

import numpy as np
from scipy.io import loadmat
from scipy.ndimage import gaussian_filter
import geo
import matplotlib.pylab as plt

isrwater = False
isnewsta = True

name = ''
lst = ''
if isrwater: name += '_rwater'
else: name += '_culled'
if isnewsta: 
    lst += '_newsta'
else: 
    lst += '_oldsta'
name += lst


def sortAB(A,B):
    B = [b for _,b in sorted(zip(A,B))]
    return B

def plotRFv(ax, data, rf_time, offsetx, norm=1):
     data[np.where(np.isnan(data))] = 0
     x1, x2 = data.copy()*norm, data.copy()*norm
     thre = 0.2*np.max(np.abs(data*norm))
     x1[np.where(data<=thre/norm)[0]] = thre
     x2[np.where(data>-thre/norm)[0]] = -thre
     ax.fill_betweenx(rf_time, thre+offsetx, x1+offsetx, facecolor='red',alpha=0.5)
     ax.fill_betweenx(rf_time, -thre+offsetx, x2+offsetx, facecolor='blue',alpha=0.5)
     ax.plot(data*norm+offsetx,rf_time,'k-',linewidth=0.3,alpha=0.3)

ccpdat = 'ccp_0.1-0.5hz%s.dat'%name
ccp = np.loadtxt(ccpdat,usecols=(2,3,4))
dist = np.unique(ccp[:,0])
ddis = np.diff(dist)[0]
dep = np.unique(ccp[:,1])[:-2]
ddep = np.diff(dep)[0]
amp = np.zeros([dep.shape[0],dist.shape[0]])

for i in range(ccp.shape[0]):
    dep_id = int(ccp[i,1]/ddep)
    dist_id = int(ccp[i,0]/ddis)
    if dep_id+1 > amp.shape[0] or dist_id+1 > amp.shape[1]: continue
    if np.isnan(ccp[i,2]): continue
    amp[dep_id,dist_id] = ccp[i,2]

amp = gaussian_filter(amp,[0.5,4])
amp[np.where(np.abs(amp)<0.1)] = 0
fig = plt.figure(figsize=(12,6))
ax = fig.add_subplot(111)
ax.pcolor(dist,dep,amp,cmap=plt.cm.seismic,vmin=-0.5,vmax=0.5)

norm = 5
line = np.array([[-133.55,-8],[-132.45,-4]])
mat = loadmat('RFdepth%s.mat'%name)['RFdepth'][0]

for i in range(len(mat)):    
    dep = mat[i]['Depthrange'][0]
    moveout = mat[i]['moveout_correct']
    projlat = mat[i]['projectlat']
    projlon = mat[i]['projectlon']
    dep, moveout, projlat, projlon = map(np.array,[dep, moveout, projlat, projlon])
    nRF = moveout.shape[1]
    gcarc,_ = geo.distance(projlat,projlon,line[0,1],line[0,0])
    dist = geo.deg2km(gcarc)
    #xoffset = dist+norm*moveout
    for j in range(nRF):
        #ax.plot(xoffset[:,j],dep,'k-',linewidth=0.2)
        plotRFv(ax,moveout[:,j],dep,dist[:,j],norm)

stalst = np.loadtxt('sta_project%s.dat'%lst,usecols=0,dtype='str').tolist()
stadis = np.loadtxt('sta_project%s.dat'%lst,usecols=1)
stadis_dic = {}
for i in range(len(stalst)):
    stadis_dic[stalst[i]] = stadis[i]

import matlab.engine
from matlab import double
eng = matlab.engine.start_matlab()
rayp = {}
watertime = {}
waterdep = {}
with open('watertime0.5hz.dat','r') as f:
    for line in f.readlines():
        line = line.replace('\n','').split()
        rayp[line[0]] = float(line[1])
        watertime[line[0]] = float(line[2])
    f.close()
for sta in rayp.keys():
    #rayp_sta, watertime_sta = map(double,[rayp[sta],watertime[sta]])
    waterdep[sta] = float(eng.time2dep(rayp[sta],watertime[sta]))
waterdep_lst = []
stadis_lst = []

ax1 = ax.twiny()
ax1.set(xticks=stadis,xticklabels=stalst,xlim=[0,460])        
ax.set(xlabel='Distance[km]',xlim=[0,460],ylabel='Depth[km]',ylim=[0,100],yticks=np.arange(0,100,10))
ax.invert_yaxis()
ax.grid(linestyle='--')

for sta in stalst:
    ax.plot(stadis_dic[sta],waterdep[sta],'kx',linewidth=2)
    waterdep_lst.append(waterdep[sta])
    stadis_lst.append(stadis_dic[sta])
waterdep_lst = sortAB(stadis_lst,waterdep_lst)
stadis_lst.sort()

ax.plot(stadis_lst,waterdep_lst,'k-',linewidth=1)
plt.savefig('CCP%s.png'%name)
plt.show()
