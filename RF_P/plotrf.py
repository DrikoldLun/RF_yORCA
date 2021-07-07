import os, obspy
import numpy as np
import matplotlib.pylab as plt

def plotRF(ax, data, rf_time, offsety, norm=1):
    y1, y2 = data.copy()*norm, data.copy()*norm
    thre = 0.03*np.max(np.abs(data)*norm)
    y1[np.where(data<=thre/norm)[0]] = thre
    y2[np.where(data>-thre/norm)[0]] = -thre
    ax.fill_between(rf_time, thre+offsety, y1+offsety, facecolor='red',alpha=0.5)
    ax.fill_between(rf_time, -thre+offsety, y2+offsety, facecolor='blue',alpha=0.5)
    ax.plot(rf_time,data*norm+offsety,'k-',linewidth=0.5,alpha=0.5)

def baz_alignment(st,sta,tb,te,figpath):
    mul = 0.05
    #sta = st[0].stats.station
    ch = st[0].stats.channel
    fs = st[0].stats.sampling_rate
    stack = np.zeros(len(st[0].data))
    baz = np.zeros(len(st))
    rf_time = np.arange(tb,te+1./fs,1./fs)
    fig = plt.figure(figsize=(8,10))

    for i in range(len(st)+1):
        if i < len(st):
            data = st[i].data
            baz[i] = st[i].stats.sac.baz
            y = mul*baz[i]
            stack += data
        else:
            stack /= len(st)
            y = mul*(np.max(np.ceil(baz*mul)/mul)+1./mul)
            data = stack
        y1, y2 = data.copy(), data.copy()
        thre = 0.02*np.max(np.abs(data))
        y1[np.where(data<=thre)[0]] = thre
        y2[np.where(data>-thre)[0]] = -thre
        plt.fill_between(rf_time, thre+y, y1+y, facecolor='red',alpha=0.5)
        plt.fill_between(rf_time, -thre+y, y2+y, facecolor='blue',alpha=0.5)
        plt.plot(rf_time,data+y,'k-',linewidth=0.5,alpha=0.5)
    ytick, ytick_label = [], []
    for i in np.arange(np.floor(mul*np.min(baz)),np.ceil(mul*np.max(baz))+1):
        ytick.append(i)
        ytick_label.append('%d'%(int(i/mul)))
    ytick.append(y)
    ytick_label.append('stack')
    plt.xlabel('t[s]'); plt.ylabel('baz[degree]'); plt.title(sta+' RF('+ch+') alignment (number: %d)'%len(st))
    plt.yticks(ytick,ytick_label)
    plt.grid()
    plt.savefig(os.path.join(figpath,'RFalignment_baz_'+sta+'_'+ch+'.png'))
    plt.show()
    plt.close(fig)

def rayp_alignment(st,sta,tb,te,figpath):
    mul = 400.
    #sta = st[0].stats.station
    ch = st[0].stats.channel
    fs = st[0].stats.sampling_rate
    stack = np.zeros(len(st[0].data))
    rayp = np.zeros(len(st))
    rf_time = np.arange(tb,te+1./fs,1./fs)
    fig = plt.figure(figsize=(8,10))

    for i in range(len(st)+1):
        if i < len(st):
            data = st[i].data
            rayp[i] = st[i].stats.sac.user4
            y = mul*rayp[i]
            stack += data
        else:
            stack /= len(st)
            y = mul*(np.max(np.ceil(rayp*mul)/mul)+1./mul)
            data = stack
        y1, y2 = data.copy(), data.copy()
        thre = 0.02*np.max(np.abs(data))
        y1[np.where(data<=thre)[0]] = thre
        y2[np.where(data>-thre)[0]] = -thre
        plt.fill_between(rf_time, thre+y, y1+y, facecolor='red',alpha=0.5)
        plt.fill_between(rf_time, -thre+y, y2+y, facecolor='blue',alpha=0.5)
        plt.plot(rf_time,data+y,'k-',linewidth=0.5,alpha=0.5)
    ytick, ytick_label = [], []
    for i in np.arange(np.floor(mul*np.min(rayp)),np.ceil(mul*np.max(rayp))+1):
        ytick.append(i)
        ytick_label.append('%.3f'%(i/mul))   
    ytick.append(y)
    ytick_label.append('stack')
    plt.xlabel('t[s]'); plt.ylabel('p[s/km]'); plt.title(sta+' RF('+ch+') alignment (number: %d)'%len(st))
    plt.yticks(ytick,ytick_label)
    plt.grid()
    plt.savefig(os.path.join(figpath,'RFalignment_rayp_'+sta+'_'+ch+'.png'))
    plt.show()
    plt.close(fig)

def classic_plot(st,sta,tb,te,figpath):
    st = sorted(st,key=lambda x:x.stats.sac.baz)
    #sta = st[0].stats.station
    ch = st[0].stats.channel
    fs = st[0].stats.sampling_rate
    stack = np.zeros(len(st[0].data))
    baz = []
    rf_time = np.arange(tb,te+1./fs,1./fs)
    for rf in st:
        baz.append(rf.stats.sac.baz)
        stack += rf.data
    stack /= len(st)
    #norm = np.abs(stack).max()*1.5
    norm = 1

    fig = plt.figure(figsize=(7,10))
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    ax1 = fig.add_axes([0.05,0.05,0.43,0.9])
    ax2 = fig.add_axes([0.52,0.05,0.43,0.9])
    for i in range(len(st)):
        plotRF(ax1,st[i].data,rf_time,i,norm)
    plotRF(ax1,stack,rf_time,len(st)+2,norm)
    ax1.set(xlim=[-1,10],ylim=[-1,len(st)+3],yticks=np.arange(0,len(st)).tolist()+[len(st)+2],yticklabels=[],\
        xlabel=r'$Time$ $after$ $P[s]$',title='%s %dRFs (%s component)'%(sta,len(st),ch))
    ax1.grid()
    ax2.scatter(baz,np.arange(len(st)),marker='o',s=10,color='none',edgecolor='blue')
    ax2.set(xlim=[0,360],ylim=[-1,len(st)+3],yticks=np.arange(0,len(st)),yticklabels=[],\
        xticks=np.arange(0,360,60),xlabel=r'$baz[deg]$',title='Backazimuth')
    ax2.grid()
    #plt.suptitle(sta+' %d RFs'%len(st))
    plt.savefig(os.path.join(figpath,'RFalignment_'+sta+'_'+ch+'.png'))
    plt.show()
    plt.close(fig)

def stanorthcompara(stalst,RFpath,orient_file=[]):
    fig = plt.figure(figsize=(3*len(stalst),10))
    stalst.sort()
    X = {}
    for sta in stalst: X[sta] = 0
    if len(orient_file) > 0:
        with open(orient_file,'r') as f:
            for line in f.readlines():
                line = line.replace('\n','').split()
                X[line[0]] = float(line[1])
            f.close()
                
    for i in range(len(stalst)):
        rffile = os.path.join(RFpath,stalst[i],'*.R')
        st = obspy.read(rffile)
        RFtime = np.linspace(st[0].stats.sac.b,st[0].stats.sac.e,st[0].stats.npts)
        st = sorted(st,key=lambda x: x.stats.sac.baz)
        ax = fig.add_subplot(1,len(stalst),i+1)
        for tr in st:
            plotRF(ax,tr.data,RFtime,(tr.stats.sac.baz-X[stalst[i]])%360,norm=2)
        ax.set(xlim=[-1,10],xlabel='t[s]',ylim=[0,360],yticks=np.arange(0,360,60),title=stalst[i])
        if i!=0: ax.set(yticklabels=[])
        ax.grid()
    plt.savefig('staN.png')
    plt.show()

if __name__ == '__main__':
    import configparser
    config = configparser.RawConfigParser()
    parafile = 'paraRF.cfg'
    config.read(parafile)
    prefix = parafile.replace(os.path.basename(parafile),'')
    out_path = os.path.join(prefix,config.get('path', 'out_path'))
    img_path = os.path.join(out_path,config.get('path', 'image_dir'))
    RF_path = os.path.join(out_path,config.get('path', 'RF_dir'))
    orient_file = os.path.join(prefix,config.get('path', 'orient_file'))
    stalst = ['CC06','CC11','WC05','EE04','WW03']
    stanorthcompara(stalst,RF_path,orient_file)
