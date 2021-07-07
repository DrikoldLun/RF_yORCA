import sys, getopt, os
def Usage():
    print('python checkrotate.py -Ssta -ConlycalC -Lbadevtlst parafile')
    sys.exit(1)

try:
    opts, args = getopt.getopt(sys.argv[1:],'S:CL')
except:
    Usage()
if opts == []:
    Usage()
onlycalC = False
badevt = False
for op, value in opts:
    if op == '-S':
        sta = value
    elif op == '-C':
        onlycalC = True
    elif op == '-L':
        badevt = True
    else:
        Usage()
parafile = sys.argv[-1]

if not os.path.isfile(parafile):
    Usage()

import numpy as np
import obspy, configparser, glob
from obspy.taup import TauPyModel
import numpy as np
import matplotlib.pylab as plt
from scipy.signal import hilbert

config = configparser.ConfigParser()
config.read(parafile)
out_path = config.get('path','out_path')
rotate_path = os.path.join(out_path,config.get('path','rotate_dir'))
RF_path = os.path.join(out_path,config.get('path','RF_dir'))
badevtlst_file = 'PCA/evtlog/'+sta+'_badRF.lst'
freqmin = config.getfloat('para','freqmin')
freqmax = config.getfloat('para','freqmax')
RF_calcwinlen = config.getfloat('para','RF_calcwinlen')
cr_dir = os.path.join('PCA','checkrotate',sta)
if not os.path.exists(cr_dir):
    os.makedirs(cr_dir)

orient = {}
with open('../result.txt','r') as f:
    for line in f.readlines():
        line = line.replace('\n','').split()
        orient[line[0]] = float(line[1])
    f.close()

'''
def plotRF(ax, data, rf_time, offsety):
    y1, y2 = data.copy(), data.copy()
    thre = 0.03*np.max(np.abs(data))
    y1[np.where(data<=thre)[0]] = thre
    y2[np.where(data>-thre)[0]] = -thre
    ax.fill_between(rf_time, thre+offsety, y1+offsety, facecolor='red',alpha=0.5)
    ax.fill_between(rf_time, -thre+offsety, y2+offsety, facecolor='blue',alpha=0.5)
    ax.plot(rf_time,data+offsety,'k-',linewidth=0.5,alpha=0.5)

def checkRF(ax1,ax2,evttime,sta):
# evttime: e.g. 201811160326
    rf_prefix = os.path.join(RF_path,sta,'_'.join([evttime,sta,'PRF']))
    rf_R = obspy.read(rf_prefix+'.R')[0].copy()
    rf_T = obspy.read(rf_prefix+'.T')[0].copy()
    fs = rf_R.stats.sampling_rate
    tb = rf_R.stats.sac.b
    te = rf_R.stats.sac.e
    #win = [-2, 5]
    #n1, n2 = [int(tmp) for tmp in [(win[0]-tb)*fs,(win[1]-tb)*fs+1]]
    RF_time = np.linspace(tb,te,num=rf_R.stats.npts,endpoint=True)
    plotRF(ax1,rf_R.data,RF_time,1)
    plotRF(ax1,rf_T.data,RF_time,-1)
    ax1.set(xlim=[-2,5],ylim=[-3,3],yticks=[-1,1],\
        yticklabels=['RF_T','RF_R'],xlabel='t[s]')

    prefix = os.path.join(rotate_path,sta,'_'.join([evttime,sta]))
    R = obspy.read(prefix+'.R')[0].copy()
    T = obspy.read(prefix+'.T')[0].copy()
    Z = obspy.read(prefix+'.Z')[0].copy()
    for trace in [R,T,Z]:
        trace.filter('bandpass',freqmin=freqmin,freqmax=freqmax,zerophase=True)
    fs = float(R.stats.sampling_rate)
    t_P = R.stats.sac.user3
    st_time = np.arange(0,R.stats.npts)/fs
    n1, n2 = [int(tmp) for tmp in [(t_P-RF_calcwinlen)*fs,(t_P+RF_calcwinlen)*fs+1]]
    norm = np.max(np.abs(R.data[n1:n2]))
    ax2.plot(st_time,2+Z.data/norm)
    ax2.plot(st_time,R.data/norm)
    ax2.plot(st_time,-2+T.data/norm)
    ax2.vlines(t_P,-3,3,color='red')
    ax2.set(xlim=[t_P-RF_calcwinlen,t_P+RF_calcwinlen],ylim=[-3,3],yticks=[-2,0,2],\
        yticklabels=['T','R','Z'],xlabel='t[s]')

def checkSKS(ax1,ax2,evttime,sta):
    t_pre, t_post = 10, 50
    prefix = os.path.join(rotate_path,sta,'_'.join([evttime,sta]))
    R = obspy.read(prefix+'.R')[0].copy()
    T = obspy.read(prefix+'.T')[0].copy()
    Z = obspy.read(prefix+'.Z')[0].copy()
    for trace in [R,T,Z]:
        trace.filter('bandpass',freqmin=0.04,freqmax=0.125,zerophase=True)
    fs = float(R.stats.sampling_rate)
    stats = R.stats.sac
    model = TauPyModel(model="iasp91")
    arrivals = model.get_travel_times(source_depth_in_km=stats.evdp,distance_in_degree=stats.gcarc,phase_list=['SKS'])
    t_SKS = arrivals[0].time
    n1, n2 = [int(tmp) for tmp in [(t_SKS-t_pre)*fs,(t_SKS+t_post)*fs+1]]
    ax1.plot(R.data[n1:n2],T.data[n1:n2])
    ax1.axis('equal')
    ax1.set(xlabel='R',ylabel='T',xticks=[],yticks=[])
    st_time = np.arange(0,R.stats.npts)/fs
    norm = np.max(np.abs(R.data[n1:n2]))
    ax2.plot(st_time,2+Z.data/norm)
    ax2.plot(st_time,R.data/norm)
    ax2.plot(st_time,-2+T.data/norm)
    ax2.vlines(t_SKS,-3,3,color='red')
    ax2.set(xlim=[t_SKS-t_pre,t_SKS+t_post],ylim=[-3,3],yticks=[-2,0,2],\
        yticklabels=['T','R','Z'],xlabel='t[s]')
'''

def checkRayleigh(ax1,ax2,evttime,sta,onlycalC=False):
    info = []
    orient_path = '../../dataquality/event/ATaCR/rdmsd_event/orient/dat_file/'+sta+'_orient.dat'
    with open(orient_path,'r') as f:
        for line in f.readlines()[1:]:
            line = line.replace('\n','').split()
            #if float(line[4]) < 0.4:
            #    continue
            info.append([float(tmp) for tmp in line[3:]])
        f.close()
    info = np.array(info)
    dev,rej,acc = [],[],[]
    for i in range(len(info)):
        if info[i,1] < 0.4:
            rej.append([info[i,1],info[i,0]])
        else:
            acc.append([info[i,1],info[i,0]])
            dev.append(info[i,0])
    dev,rej,acc = np.array(dev),np.array(rej),np.array(acc)
    ax1.scatter(rej[:,0],rej[:,1],color='red',label='$C_{zr} < 0.4$')
    ax1.scatter(acc[:,0],acc[:,1],color='blue',label='$C_{zr} \geq 0.4$')
    ax1.hlines(orient[sta],0,1,color='black',linestyle='--',label='Orientation')
    ax1.set(xlabel='$C_{zr}$',xlim=[0,1],ylabel='$Deviation[deg]$',ylim=[0,360],title='Orientation search (Sta '+sta+')\nTotally %d Mw>=6.0 events'%len(info))
    ax1.legend(loc='upper right')

    prefix = os.path.join(rotate_path,sta,'_'.join([evttime,sta]))
    R = obspy.read(prefix+'.R')[0].copy()
    T = obspy.read(prefix+'.T')[0].copy()
    Z = obspy.read(prefix+'.Z')[0].copy()
    fs = float(R.stats.sampling_rate)
    for trace in [R,T,Z]:
        trace.filter('bandpass',freqmin=0.02,freqmax=0.04,zerophase=True)
    stats = R.stats.sac
    c = 4.
    t_surf = stats.dist/c
    n1, n2 = [int(tmp) for tmp in [(t_surf-250)*fs,(t_surf+300)*fs+1]]
    if not onlycalC:
        #ax1.plot(R.data[n1:n2],Z.data[n1:n2])
        #ax1.axis('equal')
        #ax1.set(xlabel='R',ylabel='Z',xticks=[],yticks=[])
        st_time = np.arange(0,R.stats.npts)/fs
        norm = np.max(np.abs(R.data[n1:n2]))
        ax2.plot(st_time,-2+Z.data/norm,label='Z')
        ax2.plot(st_time,R.data/norm)
        ax2.plot(st_time,2+T.data/norm)
        #ax2.set(xlim=[t_surf-20,t_surf+600],ylim=[-3,3],yticks=[-2,0,2],\
        #        yticklabels=['T','R','Z'],xlabel='$t[s]$',title=' '.join(['Event %s\n'%evt_time+'Mw %.1f'%stats.mag,'Freq %.2f-%.2f Hz'%(0.02,0.04)]))
        ax2.set(xlim=[t_surf-250,t_surf+300],ylim=[-3,3],yticks=[-2,0,2],\
                yticklabels=['Z','R','T'],xlabel='$Time\ since\ Earthquake[s]$',title=' '.join(['Rotation example (Evt %s)\n'%evt_time+'Mw %.1f'%stats.mag,'Freq %.2f-%.2f Hz'%(0.02,0.04)]))
    Z_h = np.imag(hilbert(Z.data[n1:n2]))
    Szr = -np.dot(Z_h,R.data[n1:n2])
    Szz = np.dot(Z_h,Z_h)
    Srr = np.dot(R.data[n1:n2],R.data[n1:n2])
    Czr = Szr/np.sqrt(Szz*Srr)
    Czrstar = Szr/Szz

    ax2.plot(st_time,np.imag(hilbert(R.data))/norm-2,label='R_Hilbert')
    ax2.legend(loc='upper right')
    #if not onlycalC:
    #    ax1.legend(['Rayleigh wave polarization\n($C_{zr},C_{zr}^* = %.2f,%.2f$)'%(Czr,Czrstar)],loc='upper right')
    #return Czr,Czrstar

sta_dir = os.path.join(RF_path,sta)
rf_lst = glob.glob(os.path.join(sta_dir,'*.R'))
rf_lst.sort(key=lambda rf_name: obspy.read(rf_name)[0].stats.sac.baz)
baz_lst,Czr_lst,Czrstar_lst = [],[],[]
if badevt:
    baz_lst_bad,Czr_lst_bad,Czrstar_lst_bad = [],[],[]
    badevtlst = np.loadtxt(badevtlst_file).astype('int').astype('str')
ax31, ax32 = 0, 0
for rf_name in rf_lst[:1]:
    rf_name = '../RF_data/RF50_0.10to2.00/CC06/201811142121_CC06_PRF.R'
    stats = obspy.read(rf_name)[0].stats.sac
    evt_time = os.path.basename(rf_name)[:12]
    if not onlycalC:
        fig = plt.figure(figsize=(7,4))
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        #plt.suptitle(' '.join([sta,evt_time,'mag:%.1f'%stats.mag,'gcarc:%.1f'%stats.gcarc,'baz:%.1f'%stats.baz]))
        ax1 = fig.add_subplot(121)
        ax2 = fig.add_subplot(122)
    checkRayleigh(ax1,ax2,evt_time,sta,onlycalC=onlycalC)
    #Czr, Czrstar = checkRayleigh(ax1,ax2,evt_time,sta,onlycalC=onlycalC)
    #baz_lst.append(stats.baz)
    #Czr_lst.append(Czr)
    #Czrstar_lst.append(Czrstar)
    #if badevt:
    #    if evt_time in badevtlst:
    #        baz_lst_bad.append(stats.baz)
    #        Czr_lst_bad.append(Czr)
    #        Czrstar_lst_bad.append(Czrstar)
            #if not onlycalC:
            #    plt.suptitle(' '.join([sta,evt_time+'(rejected)','mag:%.1f'%stats.mag,'gcarc:%.1f'%stats.gcarc,'baz:%.1f'%stats.baz]))
    if not onlycalC:
        plt.savefig(os.path.join(cr_dir,evt_time+'_'+sta+'.png'))
        plt.show()
        #plt.show()
        plt.close('all')
'''
plt.close('all')
fig = plt.figure(figsize=(6,10))
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)
ax1.plot(baz_lst,Czr_lst,'k.')
ax2.plot(baz_lst,Czrstar_lst,'k.')
if badevt:
    ax1.plot(baz_lst_bad,Czr_lst_bad,'r.')
    ax2.plot(baz_lst_bad,Czrstar_lst_bad,'r.')
ax1.set(ylabel=r'$C_{zr}$')
ax2.set(ylim=[0,2],xlabel='baz',ylabel=r'$C_{zr}^*$')
plt.savefig(os.path.join(cr_dir,sta+'_Cvalue.png'))
#plt.show()
'''
