#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 20:59:04 2020

@author: lun
"""

import matplotlib.pyplot as plt
import os, sys, glob, re, shutil
import obspy
import numpy as np
import getopt
from matplotlib.widgets import Button
from operator import itemgetter
import initopts
import plotrf
import copy
import configparser
sys.path.append('RFscript')
import decov
from scipy.fftpack import ifft

config = configparser.RawConfigParser()

def Usage():
    print("Usage:python PCApickrf.py -Sstation -C[R/T]comp  para.cfg")
    sys.exit(1)

def get_sac():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "S:C:")
    except:
        Usage()
    if opts == []:
        Usage()
    for op, value in opts:
        if op == '-S':
            station = value
        elif op == '-C':
            comp = value
            if comp != 'R' and comp != 'T':
                print('No such component!')
                sys.exit(1)
        else:
            Usage()
    head = sys.argv[-1]
    '''
    station = 'WC03'
    comp = 'R'
    head = 'paraRF.cfg'
    '''
    
    if os.path.isfile(head):
        config.read(head)
    else:
        Usage()
    out_path = config.get('path', 'out_path')
    image_path = os.path.join(out_path,config.get('path', 'image_dir'),station)
    path = os.path.join(out_path,config.get('path', 'RF_dir'))
    time_before = config.getint('para', 'time_beforeP') * -1
    time_after = config.getint('para', 'time_afterP')
    path = os.path.join(path,station)
    filenames = []
    for filename in glob.glob(os.path.join(path, '*_PRF.'+comp)):
        filenames.append(os.path.basename(filename).split('_')[0])
    filenames.sort()
    rffiles = obspy.read(os.path.join(path, '*_PRF.'+comp))
    #staname = rffiles[0].stats.network+'.'+rffiles[0].stats.station
    return rffiles.sort(['starttime']), filenames, path, image_path, station, time_before, time_after

def sortAB(A,B):
    B = [b for _,b in sorted(zip(A,B))]
    return B

def PCA(M):
    C, sv, F = np.linalg.svd(M.T,full_matrices=False)
    meanC = np.diag(np.mean(C,0)/abs(np.mean(C,0)))
    F = np.dot(meanC,F)
    return np.power(sv,2), F.T

class plotrffig():

    fig = plt.figure(figsize=(10, 16), dpi=60)
    axPCAcull = plt.axes([0.81, 0.91, 0.07, 0.03])
    axampcull = plt.axes([0.71, 0.91, 0.07, 0.03])
    axfinish = plt.axes([0.91, 0.91, 0.07, 0.03])
    axPlot = plt.axes([0.1, 0.91, 0.07, 0.03])
    axrePlot = plt.axes([0.2, 0.91, 0.07, 0.03])
    axdelete = plt.axes([0.3, 0.91, 0.14, 0.03])
    
    bPCAcull = Button(axPCAcull, 'PCAcull')
    bampcull = Button(axampcull, 'ampcull')
    bfinish = Button(axfinish, 'Finish')
    bplot = Button(axPlot, 'Plot RFs')
    breplot = Button(axrePlot, 'rePlot')
    bdelete = Button(axdelete, 'deletebadRFs')

    def __init__(self, opts):
        self.opts = opts
        self.ax = plt.axes([0.16, 0.05, 0.5, 0.84])
        self.ax_baz = plt.axes([0.72, 0.05, 0.23, 0.84])
        self.goodrf = np.ones(opts.evt_num).astype('int')
        self.lines = [[] for i in range(opts.evt_num)]
        self.wvfillpos = [[] for i in range(opts.evt_num)]
        self.wvfillnag = [[] for i in range(opts.evt_num)]
        self.plotwave()
        # buttons
        self.bPCAcull.on_clicked(self.PCAcull)
        self.bampcull.on_clicked(self.ampcull)
        self.bfinish.on_clicked(self.finish)
        self.bplot.on_clicked(self.plot)
        self.breplot.on_clicked(self.plotwave)
        self.bdelete.on_clicked(self.deletebadRFs)
        # select waveform
        self.fig.canvas.mpl_connect('button_press_event', self.onclick)

    def onclick(self, event):
        if event.inaxes != self.ax:
            return
        click_idx = int(np.round(event.ydata))
        if click_idx > self.opts.evt_num:
            return
        if self.goodrf[click_idx-1] == 1:
            print("Selected "+os.path.basename(self.opts.filenames[click_idx-1]))
            self.goodrf[click_idx-1] = 0
            self.wvfillpos[click_idx-1].set_facecolor('gray')
            self.wvfillnag[click_idx-1].set_facecolor('gray')
        else:
            print("Canceled "+os.path.basename(self.opts.filenames[click_idx-1]))
            self.goodrf[click_idx-1] = 1
            self.wvfillpos[click_idx-1].set_facecolor('red')
            self.wvfillnag[click_idx-1].set_facecolor('blue')
        plt.draw()
        
    def finish(self, event):
        opts = self.opts
        badidx = np.where(self.goodrf == 0)[0]
        print("%d RFs are rejected" % len(badidx))
        
        evtlogdir = os.path.join('PCA/evtlog',opts.path.split('/')[-2])
        if not os.path.exists(evtlogdir):
            os.makedirs(evtlogdir)
        lines = [[], []]
        for i in range(opts.evt_num):
            lines[self.goodrf[i]].append('%s\n'%(os.path.basename(opts.filenames[i])[:12]))
        with open(os.path.join(evtlogdir,opts.staname+'_goodRF.lst'),'w+') as f:
            f.writelines(lines[1])
            f.close()
        with open(os.path.join(evtlogdir,opts.staname+'_badRF.lst'),'w+') as f:
            f.writelines(lines[0])
            f.close()
        
        stack = np.zeros(opts.RFlength)
        with open(os.path.join(evtlogdir,opts.staname+"finallist.dat"), 'w+') as fid:
            for i in range(opts.evt_num):
                evtname = os.path.basename(opts.filenames[i])
                if self.goodrf[i] == 0:
                    #os.system('rm '+os.path.join(opts.path,evtname+'*PRF.*'))
                    print("Reject PRF of "+evtname)
                    continue        
                stack += opts.rffiles[i].data
                evla = opts.rffiles[i].stats.sac.evla
                evlo = opts.rffiles[i].stats.sac.evlo
                evdp = opts.rffiles[i].stats.sac.evdp
                dist = opts.rffiles[i].stats.sac.gcarc
                baz = opts.rffiles[i].stats.sac.baz
                rayp = opts.rffiles[i].stats.sac.user4
                mag = opts.rffiles[i].stats.sac.mag
                gauss = opts.rffiles[i].stats.sac.user0
                fid.write('%s %s %6.3f %6.3f %6.3f %6.3f %6.3f %8.7f %6.3f %6.3f\n' % (evtname, 'P', evla, evlo, evdp, dist, baz, rayp, mag, gauss))
        # shutil.copy(os.path.join(opts.path, opts.staname+"finallist.dat"), os.path.join(opts.cut_path, opts.staname+"finallist.dat"))
        stack /= opts.evt_num-len(badidx)
        time_axis = np.linspace(opts.b, opts.e, opts.RFlength)
        np.savetxt(os.path.join(evtlogdir,opts.staname+'_PCA.dat'),np.vstack((time_axis,stack)),fmt='%12.7f',newline='\n')
        sys.exit(0)

    def plot(self, event):
        opts = self.opts
        st = opts.rffiles
        filenames = opts.filenames
        print('Plotting Figure of '+opts.staname)
        st_new = []
        for i in range(opts.evt_num):
            if self.goodrf[i] == 1:
                st_new.append(opts.rffiles[i])
        #plotrf.baz_alignment(st, opts.b, opts.e, opts.image_path)
        #plotrf.rayp_alignment(st, opts.b, opts.e, opts.image_path)
        plotrf.classic_plot(st_new, opts.staname, opts.b, opts.e, opts.image_path)
        print("Figure has saved into %s" % opts.image_path)
        
    def PCAcull(self, event):
        opts = self.opts
        num_dim = 1
        n1, n2 = [int(tmp) for tmp in [(opts.PCAwin[0]-opts.b)*opts.fs,(opts.PCAwin[1]-opts.b)*opts.fs+1]]
        M = []
        for i in range(opts.evt_num):
            if self.goodrf[i] == 1:
                M.append(opts.rffiles[i].data[n1:n2])
        M = np.array(M)
        eigval, eigvec = PCA(M)
        rf_temp = np.zeros(M.shape[1])
        for i in range(M.shape[0]):
            for j in range(num_dim):
                rf_temp += eigval[j]*eigvec[i,j]*M[i]
        
        auto_temp = np.sqrt(np.dot(rf_temp,rf_temp))
        k = 0
        for i in range(opts.evt_num):
            if self.goodrf[i] == 1:
                seg = opts.rffiles[i].data[n1:n2]
                auto_seg = np.sqrt(np.dot(seg,seg))
                cc = np.dot(rf_temp,seg)/(auto_temp*auto_seg)
                if cc < opts.PCAcc:
                    self.goodrf[i] = 0
                    self.wvfillpos[i].set_facecolor('gray')
                    self.wvfillnag[i].set_facecolor('gray')
                    k += 1
        plt.draw()
        print('PCAcull rejects %d rfs'%k)
        
    def ampcull(self, event):
        opts = self.opts
        time_axis = np.linspace(opts.b, opts.e, opts.RFlength)
        norm = ifft(decov.gauss_spectrum(1./opts.fs,opts.RFlength,5)).real.max()
        k = 0
        for i in range(opts.evt_num):
            if self.goodrf[i] == 1:
                Pamp = opts.rffiles[i].data[np.where((time_axis>=0)&(time_axis<0.2))].max()
                if Pamp < 0.1*norm or Pamp > norm:
                    self.goodrf[i] = 0
                    self.wvfillpos[i].set_facecolor('gray')
                    self.wvfillnag[i].set_facecolor('gray')
                    k += 1
        plt.draw()
        print('ampcull rejects %d rfs'%k)
    
    def deletebadRFs(self, event):
        opts = copy.deepcopy(self.opts)
        for i in range(self.opts.evt_num):
            if self.goodrf[i] == 0:
                evtname = os.path.basename(self.opts.filenames[i])
                print("Reject PRF of "+evtname)
                opts.rffiles.remove(self.opts.rffiles[i])
                opts.filenames.remove(self.opts.filenames[i])
                opts.baz.remove(self.opts.baz[i])
        self.opts = opts
        self.opts.evt_num = len(self.opts.baz)
        self.goodrf = np.ones(self.opts.evt_num).astype('int')
        self.opts.ylim = [0, opts.evt_num+3]
        self.plotwave()
        
    def plotwave(self, event=[]):
        ax = self.ax
        opts = self.opts
        ax.cla()
        ax.set_xlim(self.opts.xlim)
        ax.set_ylim(self.opts.ylim)
        ax.set_ylabel("Event", fontsize=20)
        ax.set_xlabel("Time after P (s)")
        ax.set_title("R component")       
        ax.grid()
        
        bound = np.zeros(opts.RFlength)
        time_axis = np.linspace(opts.b, opts.e, opts.RFlength)
        
        rffiles = [[], []]
        filenames = [[], []]
        goodrf = [[], []]
        baz = [[], []]
        
        # badrf - front, goodrf - back
        
        for i in range(opts.evt_num):
            rffiles[self.goodrf[i]].append(opts.rffiles[i])
            filenames[self.goodrf[i]].append(opts.filenames[i])
            goodrf[self.goodrf[i]].append(self.goodrf[i])
            baz[self.goodrf[i]].append(opts.baz[i])
        for tmp in [rffiles,filenames,goodrf]:
            tmp[0] = sortAB(baz[0][:],tmp[0])
            tmp[1] = sortAB(baz[1][:],tmp[1])
        baz[0].sort()
        baz[1].sort()
        opts.rffiles = rffiles[0] + rffiles[1]
        opts.filenames = filenames[0] + filenames[1]
        self.goodrf = goodrf[0] + goodrf[1]
        opts.baz = baz[0] + baz[1]
        
        stack = np.zeros(opts.RFlength)
        stalst = []
        for i in range(opts.evt_num):
            rf = opts.rffiles[i]
            #amp_axis = rf.data*opts.enf+i+1
            #amp_axis = rf.data/1.2/np.max(rf.data)+i+1
            amp_axis = 1.5*rf.data+i+1
            self.lines[i] = ax.plot(time_axis, amp_axis, color="black", linewidth=0.2)
            if self.goodrf[i] == 0:
                self.wvfillpos[i] = ax.fill_between(time_axis, amp_axis, bound+i+1, where=amp_axis >i+1, facecolor='grey', alpha=0.5)
                self.wvfillnag[i] = ax.fill_between(time_axis, amp_axis, bound+i+1, where=amp_axis <i+1, facecolor='grey', alpha=0.5)
            elif self.goodrf[i] == 1:
                stalst.append(rf.stats.station)
                self.wvfillpos[i] = ax.fill_between(time_axis, amp_axis, bound+i+1, where=amp_axis >i+1, facecolor='red', alpha=0.5)
                self.wvfillnag[i] = ax.fill_between(time_axis, amp_axis, bound+i+1, where=amp_axis <i+1, facecolor='blue', alpha=0.5)
                stack += rf.data
        
        stack = 1.5*stack/len(baz[1])+opts.evt_num+2
        ax.fill_between(time_axis, stack, bound+opts.evt_num+2, where=stack >opts.evt_num+2, facecolor='red', alpha=0.5)
        ax.fill_between(time_axis, stack, bound+opts.evt_num+2, where=stack <opts.evt_num+2, facecolor='blue', alpha=0.5)
        ax.set_yticks(np.arange(1,opts.evt_num+1).tolist()+[opts.evt_num+2])
        ax.set_yticklabels(['']*opts.evt_num+['stack'])
        self.plotbaz()
        if opts.staname != 'all':
            self.fig.suptitle("%s (Latitude: %5.2fN, Longitude: %5.2fE)" % (opts.staname, opts.stla, opts.stlo), fontsize=20)
        else:
            self.fig.suptitle("%s %d rfs %d stas" % (opts.staname, sum(self.goodrf), len(np.unique(stalst))),fontsize=20)
    def plotbaz(self):
        ax_baz = self.ax_baz
        ax_baz.cla()
        ax_baz.set_xlabel("Backazimuth")
        ax_baz.grid()
        ax_baz.scatter(self.opts.baz, np.arange(self.opts.evt_num)+1)
        ax_baz.set_xlim(0, 360)
        ax_baz.set_xticks(np.arange(0, 361, 60))
        ax_baz.set_ylim(self.opts.ylim)
        ax_baz.set_yticks(np.arange(1,self.opts.evt_num+1))
        ax_baz.set_yticklabels([])

def main():
    opts = initopts.opts()
    opts.maxidx = 20
    opts.enf = 6
    opts.rffiles, opts.filenames, opts.path, opts.image_path, opts.staname, opts.b, opts.e = get_sac()
    opts.evt_num = len(opts.rffiles)
    opts.xlim = [-2, 10]
    opts.ylim = [0, opts.evt_num+3]
    rf = opts.rffiles[0]
    opts.b = rf.stats.sac.b
    opts.e = rf.stats.sac.e
    opts.stla = rf.stats.sac.stla
    opts.stlo = rf.stats.sac.stlo
    opts.RFlength = rf.data.shape[0]
    opts.fs = float(rf.stats.sampling_rate)
    opts.PCAwin = [-1,4]
    opts.PCAcc = 0.7
    bazi = [tr.stats.sac.baz for tr in opts.rffiles]
    tmp_filenames = [[opts.filenames[i], bazi[i]] for i in range(opts.evt_num)]
    tmp_filenames = sorted(tmp_filenames, key=itemgetter(1))
    opts.filenames = [file[0] for file in tmp_filenames]
    opts.baz = np.sort(bazi)
    opts.idx_bazi = np.argsort(bazi)
    opts.rffiles = [opts.rffiles[i] for i in opts.idx_bazi]
    plotrffig(opts)

if __name__ == "__main__":
    main()
    plt.show()
