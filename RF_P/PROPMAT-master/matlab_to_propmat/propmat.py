#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 17:24:19 2020

@author: lun
"""

import matlab.engine
from matlab import double
import numpy as np
import matplotlib.pylab as plt
import matplotlib.patches as patches
from scipy.fftpack import fft, ifft
#from joblib import Parallel, delayed
import decov, sys

eng = matlab.engine.start_matlab()

def propmatsyn(zlayt,zlayb,Vp,Vs,rho,rayp,freq=[0.1,2],samprate=10.,\
               synthperiod=0.5,isplot=0,Pwin=[10.,30.],RFwin=[10.,30.],ID='example'):
    
    inc = np.rad2deg(np.arcsin(rayp*Vs[-1]))
    
    tr_time = np.linspace(-Pwin[0],Pwin[1],\
         num=int(samprate*(Pwin[1]+Pwin[0])+1),endpoint=True)
    RF_time = np.linspace(-RFwin[0],RFwin[1],\
         num=int(samprate*(RFwin[1]+RFwin[0])+1),endpoint=True)
    
    zlayt,zlayb,Vp,Vs,rho,freq,Pwin,RFwin = \
         map(double,[zlayt,zlayb,Vp,Vs,rho,freq,Pwin,RFwin])
         
    samprate,inc,synthperiod = map(float,[samprate,inc,synthperiod])

    args = ('zlayt',zlayt,'zlayb',zlayb,'Vp',Vp,'Vs',Vs,\
        'rho',rho,'samprate',samprate,'inc',inc,\
        'synthperiod',synthperiod,'freq',freq,\
        'isplot',isplot,'Pwin',Pwin,'RFwin',RFwin,'ID',ID)

    R,Z,RF = eng.propmat_syn(*args,nargout=3)
    R = np.array(R).reshape(1,-1)[0]
    Z = np.array(Z).reshape(1,-1)[0]
    RF = np.array(RF).reshape(1,-1)[0]

    return R,Z,RF,tr_time,RF_time

def noise_RF(R,Z,fs,RFwin=[10.,30.]):
    R += np.random.normal(loc=0.,scale=0.1*R.max(),size=len(R))
    Z += np.random.normal(loc=0.,scale=0.1*Z.max(),size=len(R))
    RF = decov.decovt(Z,R,1./fs,RFwin[0],RFwin[1],alpha=5)
    return RF
    
def labmodel(Mohodep,LABdep,rayp,freq=[0.1,2],Vp=[],Vs=[],rho=[],others={}):
    #Vp, Vs, rho should include 3 elements (crust,upper mantle,asthenosphere)
    #others - dictionary
    if Vp == []: Vp = [6.5, 8, 7.2]
    if Vs == []: Vs = [3.7, 4.4, 4.1]
    if rho == []: rho = [2.85, 3.4, 3.2]
    
    if Mohodep == 0:
        kwargs = {'zlayt':[0, LABdep],\
          'zlayb':[LABdep, 100],\
          'Vs':Vs[1:],\
          'Vp':Vp[1:],\
          'rho':rho[1:],\
          'rayp':rayp,\
          'freq':freq}
    elif LABdep == 0:
        kwargs = {'zlayt':[0, Mohodep],\
          'zlayb':[Mohodep, 100],\
          'Vs':Vs[:2],\
          'Vp':Vp[:2],\
          'rho':rho[:2],\
          'rayp':rayp,\
          'freq':freq}
    else:
        kwargs = {'zlayt':[0, Mohodep, LABdep],\
          'zlayb':[Mohodep, LABdep, 100],\
          'Vs':Vs,\
          'Vp':Vp,\
          'rho':rho,\
          'rayp':rayp,\
          'freq':freq}
    
    kwargs.update(others)
    return kwargs

def tf_cal(rf_cleaned,rf_original,c=0):
    # freq spectrum of RFs
    f_rf_c = fft(rf_cleaned)
    f_rf_o = fft(rf_original)
    
    # numerator & denominator of transfer function
    tf_t = f_rf_c*f_rf_o.conjugate()
    tf_b = f_rf_o*f_rf_o.conjugate()
    
    # water level
    wl = np.max(tf_b.real)*c 
    tf_b[np.where(tf_b<wl)[0]] = wl
    
    # transfer function
    tf = tf_t/tf_b
    return tf

def tf_moho(Moho_h,Vs,rayp,Vp=6.5,freq=[0.1,2],TFtype=1,others={}):
    #TFtype: 1-Moho 2-crustal multiples
    kwargs_thisfunc = {'rayp':rayp,\
                       'Vs':[Vs,4.4,4.1],\
                       'Vp':[Vp,8,7.2],\
                       'freq':freq,\
                       'others':others}
    
    kwargs_nomoho = labmodel(0,0,**kwargs_thisfunc)
    kwargs_moho = labmodel(Moho_h,0,**kwargs_thisfunc)
    
    _,_,RF_nomoho,_,RF_time = propmatsyn(**kwargs_nomoho)
    _,_,RF_moho,_,_ = propmatsyn(**kwargs_moho)
    
    if TFtype == 1:
        tf = tf_cal(RF_nomoho,RF_moho,0)
    elif TFtype == 2:
        RF_rmul = RF_moho.copy() 
        RF_rmul[np.where(RF_time>1.5)] = 0
        tf = tf_cal(RF_rmul,RF_moho,0)
    else:
        raise ValueError('TFtype could only be 1-Moho or 2-crustal multiples!')
    
    return tf

def rf_moho(Moho_h,Vs,rayp,Vp=6.5,freq=[0.1,2],others={}):
    #TFtype: 1-Moho 2-crustal multiples
    kwargs_thisfunc = {'rayp':rayp,\
                       'Vs':[Vs,4.4,4.1],\
                       'Vp':[Vp,8,7.2],\
                       'freq':freq,\
                       'others':others}
    
    kwargs_moho = labmodel(Moho_h,0,**kwargs_thisfunc)
    _,_,RF_moho,_,_ = propmatsyn(**kwargs_moho)
    return RF_moho

def test_plot(Mohodep,LABdep,rayp,freq=[0.1,2],tf_type='Moho',others={}):
    # tf_cal: Moho, Moho+LAB, Crustal multiple
    kwargs_0_0 = labmodel(0,0,rayp,freq,others=others)
    kwargs_m_0 = labmodel(Mohodep,0,rayp,freq,others=others)
    kwargs_0_l = labmodel(0,LABdep,rayp,freq,others=others)
    kwargs_m_l = labmodel(Mohodep,LABdep,rayp,freq,others=others)
    
    _,_,RF_0_0,_,RF_time = propmatsyn(**kwargs_0_0)
    _,_,RF_m_0,_,_ = propmatsyn(**kwargs_m_0)
    _,_,RF_0_l,_,_ = propmatsyn(**kwargs_0_l)
    _,_,RF_m_l,_,_ = propmatsyn(**kwargs_m_l)

    RF_0_0 *= np.abs(RF_m_l).max()/np.abs(RF_0_0).max()
    RF_m_0 *= np.abs(RF_m_l).max()/np.abs(RF_m_0).max()
    RF_0_l *= np.abs(RF_m_l).max()/np.abs(RF_0_l).max()
    
    fig = plt.figure(figsize=(15,9))
    ax1 = fig.add_subplot(231)
    ax2 = fig.add_subplot(232)
    ax3 = fig.add_subplot(233)
    ax4 = fig.add_subplot(234)
    ax5 = fig.add_subplot(235)
    ax6 = fig.add_subplot(236)
    
    if tf_type == 'Moho':
        # crustal phases removal
        tf = tf_cal(RF_0_0,RF_m_0,0)
        #tf = tf_moho(Mohodep,3.7,rayp)
        # 0/0 vs m/0 -> tranfer function
        ax1.plot(RF_time,RF_0_0,label='0/0')
        ax1.plot(RF_time,RF_m_0,label='%.1f/0'%Mohodep)
        
    elif tf_type == 'Moho+LAB':
        # tf calculated from both Moho & LAB
        kwargs_0_l0 = labmodel(0,70,rayp,freq,others=others)
        kwargs_m_l0 = labmodel(Mohodep,70,rayp,freq,others=others)
        _,_,RF_0_l0,_,_ = propmatsyn(**kwargs_0_l0)
        _,_,RF_m_l0,_,_ = propmatsyn(**kwargs_m_l0)
        tf = tf_cal(RF_0_l0,RF_m_l0,0)
        ax1.plot(RF_time,RF_0_l0,label='0/70')
        ax1.plot(RF_time,RF_m_l0,label='%.1f/70'%Mohodep)
        
    elif tf_type == 'Crustal multiple':
        # crustal multiple removal
        RF_rmul = RF_m_0.copy() 
        RF_rmul[np.where(RF_time>1.5)] = 0
        tf = tf_cal(RF_rmul,RF_m_0,0)
        ax1.plot(RF_time,RF_rmul,label='%.1f/0 (multiples) removed'%Mohodep)
        ax1.plot(RF_time,RF_m_0,label='%.1f/0'%Mohodep)
    
    RF_cleaned = ifft(fft(RF_m_l)*tf).real
    
    ax1.set(title='synthetics for calculating TF')
    
    # LAB phase & Moho effect
    ax2.plot(RF_time,RF_m_0+1,label='%.1f/0'%Mohodep)
    ax2.plot(RF_time,RF_0_l,label='0/%.1f'%LABdep)
    ax2.plot(RF_time,RF_m_l-1,label='%.1f/%.1f'%(Mohodep,LABdep))
    ax2.set(title='Moho & LAB effects',yticks=[-1,0,1],\
            yticklabels=['Moho+LAB','LAB','Moho'])
    
    # only Moho vs. overall
    ax3.plot(RF_time,RF_m_0,label='%.1f/0'%Mohodep)
    ax3.plot(RF_time,RF_m_l,label='%.1f/%.1f'%(Mohodep,LABdep))
    ax3.set(title='only Moho vs. overall')
    
    # only LAB vs. overall
    ax4.plot(RF_time,RF_0_l,label='0/%.1f'%LABdep)
    ax4.plot(RF_time,RF_m_l,label='%.1f/%.1f'%(Mohodep,LABdep))
    ax4.set(title='only LAB vs. overall')
    
    # overall vs. cleaned RF
    ax5.plot(RF_time,RF_m_l,label='%.1f/%.1f'%(Mohodep,LABdep))
    ax5.plot(RF_time,RF_cleaned,label='Moho-cleaned')
    ax5.set(title='overall vs. cleaned RF')
    
    # only LAB vs. cleaned RF
    ax6.plot(RF_time,RF_0_l,label='0/%.1f'%LABdep)
    ax6.plot(RF_time,RF_cleaned,label='Moho-cleaned')
    ax6.set(title='only LAB vs. cleaned RF')
    
    for ax in [ax1,ax2,ax3,ax4,ax5,ax6]:
        ax.set(xlim=[-2,10],xlabel='t[s]')
        ax.legend(loc='upper right')
        
    plt.suptitle('Moho:%.1fkm LAB:%.1fkm\nFreq:%.1f-%.1fHz\ntf_type:%s'%(Mohodep,LABdep,freq[0],freq[1],tf_type))
    
    #plt.savefig('/home/lun/Desktop/syn_test_crustclean/1hz/'+'%.1f_%.1f_%s.png'%(Mohodep,LABdep,tf_type))
    plt.show()

def plotRF(ax, data, rf_time, offsety, norm=1):
    y1, y2 = data.copy()/norm, data.copy()/norm
    thre = 0.03*np.max(np.abs(data))
    y1[np.where(data<=thre)[0]] = thre
    y2[np.where(data>-thre)[0]] = -thre
    ax.fill_between(rf_time, thre+offsety, y1+offsety, facecolor='red',alpha=0.5)
    ax.fill_between(rf_time, -thre+offsety, y2+offsety, facecolor='blue',alpha=0.5)
    ax.plot(rf_time,data/norm+offsety,'k-',linewidth=0.5,alpha=0.5)

def plot_rfsclean(H,Vs,rayp=0,sta='',samprate=10,parafile='../../paraRF.cfg',ID='example',evtculling=False,TFtype=2,isshow=False):
    import obspy, os, configparser, glob
    if rayp == 0:
        uniformrayp = False
        file_suffix = ''
    else:
        uniformrayp = True
        file_suffix = '(rayp=%.3f)'%rayp

    # path and input trs
    config = configparser.RawConfigParser()
    config.read(parafile)
    prefix = parafile.replace(os.path.basename(parafile),'')
    out_path = os.path.join(prefix,config.get('path', 'out_path'))
    img_path = os.path.join(out_path,config.get('path', 'image_dir'),sta)
    RF_path = os.path.join(out_path,config.get('path', 'RF_dir'),sta)
    CCPdat_path = os.path.join(out_path,'CCP',config.get('path', 'RF_dir'),'RFdat_rmul',sta)
    if not os.path.exists(CCPdat_path): os.makedirs(CCPdat_path)
    evtlst_path = os.path.join(prefix,config.get('path', 'cull_evtlist'),config.get('path', 'RF_dir'))
    finallist = os.path.join(evtlst_path,'XX.%sfinallist.dat'%sta)
    os.system('cp %s %s'%(finallist,CCPdat_path))

    st = []
    rffilenames = []
    RF_names = glob.glob(os.path.join(RF_path,'*.R'))
    if evtculling:
        evtlist_file = os.path.join(evtlst_path,'XX.%s_goodRF.lst'%sta)
        evtlist = np.loadtxt(evtlist_file).astype('int').astype('str')
        for RF_name in RF_names:
            evt_time = os.path.basename(RF_name)[:12]
            if evt_time in evtlist:
                st.append(obspy.read(RF_name)[0].copy())
                rffilenames.append(os.path.basename(RF_name))
        file_suffix += '_culled'
        ID += '_culled'
    else:
        for RF_name in RF_names:
            st.append(obspy.read(RF_name)[0].copy())
            rffilenames.append(os.path.basename(RF_name))
    file_suffix += '%d'%TFtype

    baz = [tr.stats.sac.baz for tr in st]
    id_sort = np.argsort(baz)
    rffilenames = [rffilenames[i] for i in id_sort]
    baz = [baz[i] for i in id_sort]

    st = sorted(st,key=lambda x:x.stats.sac.baz)
    
    # TF clean
    RF_nomoho = rf_moho(0,0,0.06,others={'ID':ID})
    tb = st[0].stats.sac.b
    te = st[0].stats.sac.e
    fs = st[0].stats.sampling_rate
    RF_time = np.linspace(tb,te,num=int((te-tb)*samprate)+1,endpoint=True)

    st_clean = []
    # RF_rmul, numerator of TF
    if TFtype == 1: RF_rmul = RF_nomoho #Moho

    if uniformrayp:
        RF_moho = rf_moho(H,Vs,rayp,others={'ID':ID})
        if TFtype == 2:
            RF_rmul = RF_moho.copy()
            RF_rmul[np.where(RF_time>1.5)] = 0
        for i in range(len(st)):
            rf = st[i]
            Pamp = rf.data[np.where((RF_time>-0.2)&(RF_time<0.2))].max()
            fs = rf.stats.sampling_rate
            rf.decimate(factor=int(fs//samprate),no_filter=True,strict_length=False)
            tf = tf_cal(RF_rmul/np.abs(RF_rmul).max()*Pamp,RF_moho/np.abs(RF_moho).max()*Pamp)
            rf_clean = ifft(fft(rf.data)*tf).real
            st_clean.append(rf_clean)
            np.savetxt(os.path.join(CCPdat_path,rffilenames[i]),rf_clean,newline='\n')
    else:
        for i in range(len(st)):
            rf = st[i]
            Pamp = rf.data[np.where((RF_time>-0.2)&(RF_time<0.2))].max()
            fs = rf.stats.sampling_rate
            rf.decimate(factor=int(fs//samprate),no_filter=True,strict_length=False)
            rayp = rf.stats.sac.user4
            RF_moho = rf_moho(H,Vs,rayp,others={'ID':ID})
            if TFtype == 2:
                RF_rmul = RF_moho.copy()
                RF_rmul[np.where(RF_time>1.5)] = 0
            tf = tf_cal(RF_rmul/np.abs(RF_rmul).max()*Pamp,RF_moho/np.abs(RF_moho).max()*Pamp)
            rf_clean = ifft(fft(rf.data)*tf).real
            st_clean.append(rf_clean)
            np.savetxt(os.path.join(CCPdat_path,rffilenames[i]),rf_clean,newline='\n')
    
    fig = plt.figure(figsize=(7,10))
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    ax1 = fig.add_axes([0.05,0.05,0.43,0.9])
    ax2 = fig.add_axes([0.52,0.05,0.43,0.9])
    stack, stack_clean = np.zeros(len(RF_time)), np.zeros(len(RF_time))
    norm = 1
    for i in range(len(st)):
        ax1.plot(RF_time,norm*st[i].data+i,'b-')
        ax1.plot(RF_time,norm*st_clean[i]+i,'r-')
        stack += st[i].data
        stack_clean += st_clean[i]

    stack /= len(st)
    stack_clean /= len(st)
    ax1.plot(RF_time,norm*stack+len(st)+2,'b-')
    ax1.plot(RF_time,norm*stack_clean+len(st)+2,'r-')
    ax1.set(xlim=[-1,10],ylim=[-1,len(st)+3],yticks=np.arange(0,len(st)).tolist()+[len(st)+2],yticklabels=[],\
        xlabel=r'$Time$ $after$ $P[s]$',title='%s %dRFs'%(sta,len(st)))
    ax1.grid()

    ax2.scatter(baz,np.arange(len(st)),marker='o',s=10,color='none',edgecolor='blue')
    ax2.set(xlim=[0,360],ylim=[-1,len(st)+3],yticks=np.arange(0,len(st)),yticklabels=[],\
        xticks=np.arange(0,360,60),xlabel=r'$baz[deg]$',title='Backazimuth')
    ax2.grid()
    plt.savefig(os.path.join(img_path,'RFclean_alignment_'+sta+file_suffix+'.png'))
    if isshow: plt.show()

def HVs_singletrace(data,rayp,sta='',ID='example',Vp=6.5,samprate=10,Hlim=[4,8],Vslim=[3.2,4.2],\
                   Hsearch=[4.5,7.5],Vssearch=[3.2,4.2],TFtype=1,freq=[0.1,2],RFwin=[10,30]):
    dH, dVs = 0.05, 0.05
    nH = int((Hlim[1]-Hlim[0])/dH+1)
    nVs = int((Vslim[1]-Vslim[0])/dVs+1)
    Hbg = Hlim[0]+np.arange(0,nH)*dH
    Vsbg = Vslim[0]+np.arange(0,nVs)*dVs
    searchHb = int((Hsearch[0]-Hbg[0])/dH)
    searchHe = int((Hsearch[1]-Hbg[0])/dH)
    searchVsb = int((Vssearch[0]-Vsbg[0])/dVs)+1
    searchVse = int((Vssearch[1]-Vsbg[0])/dVs)+1
    s = np.zeros([len(Vsbg),len(Hbg)])
    RF_time = np.linspace(-RFwin[0],RFwin[1],num=len(data),endpoint=True)

    Pamp = data[np.where((RF_time>-0.2)&(RF_time<0.2))].max()
    if Pamp <=0: raise ValueError('P arrival not found')

    RF_nomoho = rf_moho(0,0,rayp,others={'ID':ID})
    
    anomaly = []
    
    for Vsid in range(len(Vsbg)):
        Vs = Vsbg[Vsid]
        for hid in range(len(Hbg)):
            h = Hbg[hid]
            print(h,Vs)
            try:
                RF_moho = rf_moho(h,Vs,rayp,others={'ID':ID})
                # normalization
                tf = tf_cal(RF_nomoho/np.abs(RF_nomoho).max()*Pamp,RF_moho/np.abs(RF_moho).max()*Pamp)
                rf_clean = ifft(fft(data)*tf).real
                rf_clean_seg = rf_clean[np.where((RF_time>0.5)&(RF_time<4))]
                s[Vsid,hid] += np.linalg.norm(rf_clean_seg)
            except:
                print('synthetic failed')
                anomaly.append([h,Vs])

    s[np.where(s==0)] = np.max(s)
    np.save('s_matrice/tmp_%s.npy'%sta,s)
    
    #s = np.load('s_matrice/tmp_%s.npy'%sta)
    #anomaly = [[5.75,3.9],[5.8,3.9],[5.85,3.9],[5.75,4.15],[5.8,4.15],[5.85,4.15]]
    '''
    # parallel
    kh_id = [(kid,hid) for kid in range(len(kbg)) for hid in range(len(Hbg))]
    def gs_hk(kh_pair):
        kid, hid = kh_pair
        k, h = kbg[kid], Hbg[hid]
        print(k,h)
        Vs = Vp/k
        tf = tf_moho(h,Vs,rayp,Vp=Vp,freq=freq,TFtype=TFtype,others={'RFwin':RFwin,'samprate':samprate})
        rf_clean = ifft(fft(data)*tf).real
        return (rf_clean,tf)
    rf_tf = Parallel(n_jobs=1)(delayed(gs_hk)(kh_pair) for kh_pair in kh_id)
    s[tuple(kh_id)] = [np.sqrt(np.linalg.norm(tmp[0])/len(data)) for tmp in rf_tf]
    '''
    # parallel
    
    s_search = s[searchVsb:searchVse,searchHb:searchHe]
    min_s = np.min(s_search)
    i_search, j_search = np.where(s_search==min_s)
    i = i_search[0]+searchVsb
    j = j_search[0]+searchHb
                
    Vsfinal = Vsbg[i]
    Hfinal = Hbg[j]
    
    RF_moho = rf_moho(Hfinal,Vsfinal,rayp,others={'ID':ID})
    # clean all crustal phases using best (H,Vs) pair
    tf = tf_cal(RF_nomoho/np.abs(RF_nomoho).max()*Pamp,RF_moho/np.abs(RF_moho).max()*Pamp)
    rf_object = ifft(fft(data)*tf).real

    # clean multiples using best (H,Vs) pair
    RF_rmul = RF_moho.copy()
    RF_rmul[np.where(RF_time>1.5)] = 0
    tf = tf_cal(RF_rmul/np.abs(RF_rmul).max()*Pamp,RF_moho/np.abs(RF_moho).max()*Pamp)
    #tf = propmat.tf_cal(RF_nomoho,RF_moho)
    rf_clean = ifft(fft(data)*tf).real
    
    print('H=%.2f[km], Vs=%.2f[km/s]'%(Hfinal,Vsfinal))
    with open('HVs.dat','a+') as f:
        f.write('%s %.2f %.2f\n'%(sta,Hfinal,Vsfinal))
        f.close()
    
    fig = plt.figure(figsize=(6,8))
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    ax1 = fig.add_axes([0.1,0.06,0.76,0.45])
    ax2 = fig.add_axes([0.1,0.57,0.8,0.17])
    ax3 = fig.add_axes([0.1,0.76,0.8,0.17])
    ax4 = fig.add_axes([0.88,0.06,0.02,0.45]) # colorbar
    
    cb = ax1.pcolor((Hbg-dH/2.).tolist()+[np.max(Hbg)+dH/2.],(Vsbg-dVs/2).tolist()+[np.max(Vsbg)+dVs/2.],\
        s,cmap=plt.cm.gist_rainbow,vmin=np.min(s)/2,vmax=2)
    ax1.add_patch(patches.Rectangle((Hsearch[0],Vssearch[0]),Hsearch[-1]-Hsearch[0],\
        Vssearch[-1]-Vssearch[0],fill=False,linewidth=1,linestyle='dashed'))
    ax1.vlines(Hfinal,Vsbg[0],Vsbg[-1],color='r',linestyle='dashed')
    ax1.hlines(Vsfinal,Hbg[0],Hbg[-1],color='r',linestyle='dashed')
    if len(anomaly) > 0:
        for pts in anomaly:
            ax1.add_patch(patches.Rectangle((pts[0]-dH/2.,pts[1]-dVs/2.),dH,dVs,fill=True,color='white',linewidth=0))
    fig.colorbar(cb,cax=ax4,label='RF norm')
    ax1.set(xlabel='H[km]',ylabel='Vs[km/s]',xlim=Hlim,ylim=Vslim)

    ax2.plot(RF_time,data,'b-',label='RF_obs')
    ax2.plot(RF_time,rf_clean,'r-',label='TFtype=2')
    ax2.set(xlim=[-1,10],xlabel='t[s]')
    ax2.legend(loc='upper right')

    ax3.plot(RF_time,data,'b-',label='RF_obs')
    ax3.plot(RF_time,rf_object,'r-',label='TFtype=1')
    ax3.set(title='%s avg_rayp:%.3f'%(sta,rayp),xlim=[-1,10],xticklabels=[])
    ax3.legend(loc='upper right')

    plt.savefig('HVs'+sta+'.png')
    plt.show()

def gs_station(sta,rayp=0,ID='example',Vp=6.5,samprate=10,Hlim=[4,8],Vslim=[3.,4.2],Hsearch=[4.5,7.5],\
               Vssearch=[3,4.2],freq=[0.1,2],parafile='../../paraRF.cfg',evtculling=False):
    import obspy, os, configparser, glob
    if rayp == 0:
        uniformrayp = False
        file_suffix = ''
    else:
        uniformrayp = True
        file_suffix = '(rayp=%.3f)'%rayp
    config = configparser.RawConfigParser()
    config.read(parafile)
    prefix = parafile.replace(os.path.basename(parafile),'')
    out_path = os.path.join(prefix,config.get('path', 'out_path'))
    img_path = os.path.join(out_path,config.get('path', 'image_dir'),sta)
    RF_path = os.path.join(out_path,config.get('path', 'RF_dir'),sta)
    
    dH, dVs = 0.1, 0.1
    nH = int((Hlim[1]-Hlim[0])/dH+1)
    nVs = int((Vslim[1]-Vslim[0])/dVs+1)
    Hbg = Hlim[0]+np.arange(0,nH)*dH
    Vsbg = Vslim[0]+np.arange(0,nVs)*dVs
    #Hbg = np.linspace(Hlim[0],Hlim[1],nH,endpoint=True)
    #Vsbg = np.linspace(Vslim[0],Vslim[1],nVs,endpoint=True)
    searchHb = int((Hsearch[0]-Hbg[0])/dH)
    searchHe = int((Hsearch[1]-Hbg[0])/dH)
    searchVsb = int((Vssearch[0]-Vsbg[0])/dVs)+1
    searchVse = int((Vssearch[1]-Vsbg[0])/dVs)+1
    s = np.zeros([len(Vsbg),len(Hbg)])
    
    if evtculling:
        evtlist_file = os.path.join(prefix,config.get('path', 'cull_evtlist',config.get('path', 'RF_dir')),'%s_goodRF.lst'%sta)
        evtlist = np.loadtxt(evtlist_file).astype('int').astype('str')
        RF_names = glob.glob(os.path.join(RF_path,'*.R'))
        st = []
        for RF_name in RF_names:
            evt_time = os.path.basename(RF_name)[:12]
            if evt_time in evtlist:
                st.append(obspy.read(RF_name)[0])
        file_suffix += '_culled'
        ID += '_culled'
    else:
        st = obspy.read(os.path.join(RF_path,'*_PRF.R'))
    RF_nomoho = rf_moho(0,0,0.06,others={'ID':ID})
    tb = st[0].stats.sac.b
    te = st[0].stats.sac.e
    fs = st[0].stats.sampling_rate
    RF_time = np.linspace(tb,te,num=int((te-tb)*samprate)+1,endpoint=True)
    
    for Vsid in range(len(Vsbg)):
        Vs = Vsbg[Vsid]
        for hid in range(len(Hbg)):
            h = Hbg[hid]
            print(h,Vs)
            i = 0
            if uniformrayp:
                try:
                    RF_moho = rf_moho(h,Vs,rayp,others={'ID':ID})
                except:
                    print('synthetic failed, continue')
                    continue
                for rf in st:
                    rf = rf.copy()
                    rf.decimate(factor=int(fs//samprate),no_filter=True,strict_length=False)
                    Pamp = rf.data[np.where((RF_time>-0.2)&(RF_time<0.2))].max()
                    if Pamp <= 0: 
                        print('P arrival not found, continue')
                        continue
                    # normalization
                    tf = tf_cal(RF_nomoho/np.abs(RF_nomoho).max()*Pamp,RF_moho/np.abs(RF_moho).max()*Pamp)
                    #tf = tf_moho(h,Vs,rayp,Vp=Vp,freq=freq,TFtype=TFtype,others={'RFwin':[-tb,te],'samprate':samprate})
                    rf_clean = ifft(fft(rf.data)*tf).real
                    rf_clean_seg = rf_clean[np.where((RF_time>0.5)&(RF_time<4))]
                    s[Vsid,hid] += np.linalg.norm(rf_clean_seg)
                    i += 1
            else:
                for rf in st:
                    rf = rf.copy()
                    rf.decimate(factor=int(fs//samprate),no_filter=True,strict_length=False)
                    Pamp = rf.data[np.where((RF_time>-0.2)&(RF_time<0.2))].max()
                    if Pamp <=0: 
                        print('P arrival not found, continue')
                        continue
                    rayp = rf.stats.sac.user4
                    
                    '''
                    RF_moho = rf_moho(h,Vs,rayp,others={'ID':ID})
                    # normalization
                    tf = tf_cal(RF_nomoho/np.abs(RF_nomoho).max()*Pamp,RF_moho/np.abs(RF_moho).max()*Pamp)
                    #tf = tf_moho(h,Vs,rayp,Vp=Vp,freq=freq,TFtype=TFtype,others={'RFwin':[-tb,te],'samprate':samprate})
                    rf_clean = ifft(fft(rf.data)*tf).real
                    rf_clean_seg = rf_clean[np.where((RF_time>0.5)&(RF_time<4))]
                    s[Vsid,hid] += np.linalg.norm(rf_clean_seg)
                    '''
                    
                    try:
                        RF_moho = rf_moho(h,Vs,rayp,others={'ID':ID})
                        # normalization
                        tf = tf_cal(RF_nomoho/np.abs(RF_nomoho).max()*Pamp,RF_moho/np.abs(RF_moho).max()*Pamp)
                        #tf = tf_moho(h,Vs,rayp,Vp=Vp,freq=freq,TFtype=TFtype,others={'RFwin':[-tb,te],'samprate':samprate})
                        rf_clean = ifft(fft(rf.data)*tf).real
                        rf_clean_seg = rf_clean[np.where((RF_time>0.5)&(RF_time<4))]
                        s[Vsid,hid] += np.linalg.norm(rf_clean_seg)
                        i += 1
                    except:
                        print('synthetic failed, continue')
                        continue
                    
            if i==0: continue
            s[Vsid,hid] /= i
                
    s[np.where(s==0)] = np.max(s)
    np.save(os.path.join(img_path,'s'+file_suffix+'.npy'),s)
    
    s_search = s[searchVsb:searchVse,searchHb:searchHe]
    min_s = np.min(s_search)
    i_search, j_search = np.where(s_search==min_s)
    i = i_search[0]+searchVsb
    j = j_search[0]+searchHb
                
    Vsfinal = Vsbg[i]
    Hfinal = Hbg[j]

    print('H=%.2f[km], Vs=%.2f[km/s]'%(Hfinal,Vsfinal))
    
    rf_lst = []
    rf_clean_lst = []
    if uniformrayp:
        RF_moho = rf_moho(Hfinal,Vsfinal,rayp,others={'ID':ID})
    for rf in st:
        rf = rf.copy()
        fs = rf.stats.sampling_rate
        rf.decimate(factor=int(fs//samprate),no_filter=True,strict_length=False)
        Pamp = rf.data[np.where((RF_time>-0.2)&(RF_time<0.2))].max()
        if Pamp <=0: 
            print('P arrival not found, continue')
            continue
        if not uniformrayp:
            rayp = rf.stats.sac.user4
            RF_moho = rf_moho(Hfinal,Vsfinal,rayp,others={'ID':ID})
            st.stats.channel = new_ch[i]
        RF_rmul = RF_moho.copy()
        RF_rmul[np.where(RF_time>1.5)] = 0
        tf = tf_cal(RF_rmul/np.abs(RF_rmul).max()*Pamp,RF_moho/np.abs(RF_moho).max()*Pamp)
        #tf = propmat.tf_cal(RF_nomoho,RF_moho)
        #tf = tf_moho(Hfinal,Vs,rayp,Vp=Vp,freq=freq,TFtype=2,others={'RFwin':[-tb,te],'samprate':samprate})
        rf_clean = ifft(fft(rf.data)*tf).real
        rf_lst.append(rf.data)
        rf_clean_lst.append(rf_clean)
    stack = np.array(rf_lst).mean(axis=0)
    stack_clean = np.array(rf_clean_lst).mean(axis=0)
        
    fig = plt.figure(figsize=(6,8))
    ax1 = fig.add_axes([0.05,0.05,0.85,0.6])
    ax2 = fig.add_axes([0.05,0.72,0.85,0.25])
    ax3 = fig.add_axes([0.92,0.05,0.02,0.6]) # colorbar
    
    cb = ax1.pcolor(Hbg,Vsbg,s,cmap=plt.cm.rainbow)
    ax1.add_patch(patches.Rectangle((Hsearch[0],Vssearch[0]),Hsearch[-1]-Hsearch[0],\
        Vssearch[-1]-Vssearch[0],fill=False,linewidth=1,linestyle='dashed'))
    ax1.vlines(Hfinal,Vsbg[0],Vsbg[-1],color='r',linestyle='dashed')
    ax1.hlines(Vsfinal,Hbg[0],Hbg[-1],color='r',linestyle='dashed')
    
    fig.colorbar(cb,cax=ax3,label='RF norm')
    ax1.set(title='',xlabel='H[km]',ylabel='Vs[km/s]',xlim=Hlim,ylim=Vslim)
    ax2.plot(RF_time,stack,label='original')
    ax2.plot(RF_time,stack_clean,label='cleaned')
    ax2.set(xlim=[-1,10],xlabel='t[s]')
    ax2.legend(loc='upper right')
    if uniformrayp: ax2.set(title='%.3f'%rayp)
    plt.savefig(os.path.join(img_path,'HVs'+sta+file_suffix+'.png'))
    #plt.show()
    
    with open(os.path.join(img_path,sta+file_suffix+'_HVs.txt'),'w+') as f:
        f.write('[HVspara]\n')
        f.write('station = %s\n'%sta)
        f.write('Vp = %.2f\n'%Vp)
        f.write('samprate = %.2f\n'%samprate)
        f.write('Hlimit[km] = %.1f-%.1f\n'%(Hlim[0],Hlim[1]))
        f.write('Vslimit = %.2f-%.2f\n'%(Vslim[0],Vslim[1]))
        f.write('H_search_range[km] = %.1f-%.1f\n'%(Hsearch[0],Hsearch[1]))
        f.write('Vs_search_range = %.2f-%.2f\n'%(Vssearch[0],Vssearch[1]))
        #f.write('TFtype = %d\n'%TFtype)
        f.write('freq = %.2f,%.2f\n'%(freq[0],freq[1]))
        f.write('[Result]\n')
        f.write('H = %.2f\n'%Hfinal) 
        f.write('Vs = %.2f\n'%Vsfinal)
        f.close()

def Usage():
    print('python propmat.py -Ssta')
    sys.exit(1)

if __name__ == '__main__':
    
    import getopt, configparser, os
    try:
        opts, args = getopt.getopt(sys.argv[1:],'S:')
    except:
        Usage()
    if opts == []: Usage()

    for op, value in opts:
        if op == '-S':
            sta = value
        else:
            Usage()
    
    config = configparser.RawConfigParser()
    parafile = '../../paraRF.cfg'
    config.read(parafile)
    prefix = parafile.replace(os.path.basename(parafile),'')
    out_path = os.path.join(prefix,config.get('path', 'out_path'))
    img_path = os.path.join(out_path,config.get('path', 'image_dir'))
    RF_path = os.path.join(out_path,config.get('path', 'RF_dir'))
    evtlst_path = os.path.join(prefix,config.get('path', 'cull_evtlist'),config.get('path', 'RF_dir'))
    finallist = os.path.join(evtlst_path,'XX.%sfinallist.dat'%sta)
    #finallist = os.path.join(RF_path,sta,'XX.%sfinallist.dat'%sta)
    rayp_lst = np.loadtxt(finallist,usecols=7)
    rayp = rayp_lst.mean()
    #print(rayp)
    obs_file = os.path.join(evtlst_path,'XX.%s_PCA.dat'%sta)
    data = np.loadtxt(obs_file)[1][::5]

    HVs_singletrace(data,rayp,sta=sta,ID=sta,Vssearch=[3.5,4.0]) #use this line to get optimal crustal parameters (H and Vs) for crustal multiple cleanning


    # use the following code to correct each individual RFs
    with open('HVs.dat','r') as f:
        for line in f.readlines():
            line = line.replace('\n','').split()
            if line[0] == sta:
                Hfinal = float(line[1])
                Vsfinal = float(line[2])
                break
        f.close()
    plot_rfsclean(Hfinal,Vsfinal,rayp=rayp,sta=sta,ID=sta,evtculling=True,isshow=True)

    # no sense code
    #gs_station(sta,rayp,ID=sta)
    #gs_station(sta,ID=sta)
    #gs_station(sta,rayp,ID=sta,evtculling=True)
    #gs_station(sta,ID=sta,evtculling=True)
    '''
    rayp = np.sin(np.deg2rad(15))/7.2
    test_plot(6,50,rayp,freq=[0.1,0.5],tf_type='Moho',others={'synthperiod':4})
    '''
#rayp = np.sin(np.deg2rad(15))/4.1
#for moho in [5,5.5,6]:
#    for lab in [50,70,90]:
#        for tf_type in ['Moho','Moho+LAB','Crustal multiple']:
#            test_plot(moho,lab,rayp,freq=[0.1,1],tf_type=tf_type)
