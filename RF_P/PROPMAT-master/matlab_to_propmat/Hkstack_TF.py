import numpy as np
import configparser, sys, os, obspy, tkinter
from joblib import Parallel, delayed
import matplotlib.pylab as plt
import matplotlib.patches as patches
import propmat
from propmat import tf_moho
from scipy.fftpack import fft, ifft

config = configparser.RawConfigParser()

def elipab(o,lx,ly):
    theta = np.linspace(0,2*np.pi,20,endpoint=True)
    x = lx*np.cos(theta)+o[0]
    y = ly*np.sin(theta)+o[1]
    return x,y

def Hkstack(sta,Vp,samprate,Hlim,klim,Hsearch,ksearch,TFtype,freq,RF_path,img_path):
    dH, dk = 0.1, 0.01
    nH = (Hlim[1]-Hlim[0])/dH+1
    nk = (klim[1]-klim[0])/dk+1
    Hbg = np.linspace(Hlim[0],Hlim[1],nH,endpoint=True)
    kbg = np.linspace(klim[0],klim[1],nk,endpoint=True)
    searchHb = int((Hsearch[0]-Hbg[0])/dH)
    searchHe = int((Hsearch[1]-Hbg[0])/dH)
    searchkb = int((ksearch[0]-kbg[0])/dk)+1
    searchke = int((ksearch[1]-kbg[0])/dk)+1
    s = np.zeros([len(kbg),len(Hbg)])

    print(1)
    path = os.path.join(RF_path,sta)
    st = obspy.read(os.path.join(path,'*_PRF.R'))[:1]
    for kid in range(len(kbg)):
        k = kbg[kid]
        for hid in range(len(Hbg)):
            h = Hbg[hid]
            i = 0
                
            for rf in st:
                rf = rf.copy()
                fs = rf.stats.sampling_rate
                rf.decimate(factor=int(fs//samprate),no_filter=True,strict_length=False)
                tb = rf.stats.sac.b
                te = rf.stats.sac.e
                rayp = rf.stats.sac.user4
                Vs = Vp/k
                tf = propmat.tf_moho(h,Vs,rayp,Vp=Vp,freq=freq,TFtype=TFtype,others={'RFwin':[-tb,te],'samprate':samprate})
                rf_clean = ifft(fft(rf.data)*tf).real
                s[kid,hid] += np.sqrt(np.dot(rf_clean-rf.data,rf_clean-rf.data)/len(rf_clean))
                i += 1
                print(i)
            
    print(2)
    s /= len(st)

    s_search = s[searchkb:searchke,searchHb:searchHe]
    max_s = np.max(s_search)
    i_search, j_search = np.where(s_search==max_s)
    i = i_search[0]+searchkb
    j = j_search[0]+searchHb
                
    kfinal = kbg[i]
    Hfinal = Hbg[j]

    print('H=%.2f[km], k=%.2f'%(Hfinal,kfinal))

    sigs = np.std(s.reshape(1,-1)[0])
    pdk2 = np.diff(s[i-1:i+2,j],2)/np.power(dk,2)
    pdh2 = np.diff(s[i,j-1:j+2],2)/np.power(dH,2)
    sigk = np.sqrt(2*sigs/np.abs(pdk2))
    sigh = np.sqrt(2*sigs/np.abs(pdh2))
    
    rf_lst = []
    rf_clean_lst = []
    for rf in st:
        rf = rf.copy()
        fs = rf.stats.sampling_rate
        rf.decimate(factor=int(fs//samprate),no_filter=True,strict_length=False)
        tb = rf.stats.sac.b
        te = rf.stats.sac.e
        rayp = rf.stats.sac.user4
        Vs = Vp/kfinal
        tf = propmat.tf_moho(Hfinal,Vs,rayp,Vp=Vp,freq=freq,TFtype=TFtype,others={'RFwin':[-tb,te],'samprate':samprate})
        rf_clean = ifft(fft(rf.data)*tf).real
        rf_lst.append(rf.data)
        rf_clean_lst.append(rf_clean)
    stack = np.array(rf_lst).mean(axis=0)
    stack_clean = np.array(rf_clean_lst).mean(axis=0)
        
    fig = plt.figure(figsize=(6,10))
    ax1 = fig.add_axes([0.05,0.05,0.9,0.6])
    ax2 = fig.add_axes([0.05,0.7,0.9,0.25])
    
    cb = ax1.pcolor(Hbg,kbg,s,cmap=plt.cm.Blues)
    ax1.add_patch(patches.Rectangle((Hsearch[0],ksearch[0]),Hsearch[-1]-Hsearch[0],\
        ksearch[-1]-ksearch[0],fill=False,linewidth=1,linestyle='dashed'))
    ax1.vlines(Hfinal,kbg[0],kbg[-1],color='r',linestyle='dashed')
    ax1.hlines(kfinal,Hbg[0],Hbg[-1],color='r',linestyle='dashed')
    try:
        ex, ey = elipab([Hfinal,kfinal],sigh,sigk)
        ax1.plot(ex,ey,'k')
        print('Hstd=%.2f[km], kstd=%.2f'%(sigh[0],sigk[0]))
    except:
        pass
    fig.colorbar(cb,ax=ax1)
    ax1.set(title='',xlabel='H[km]',ylabel='k',xlim=Hlim,ylim=klim)
    
    RF_time = np.linspace(tb,te,num=len(stack),endpoint=True)
    ax2.plot(RF_time,stack,label='original stack')
    ax2.plot(RF_time,stack_clean,label='cleaned stack')
    ax2.set(xlim=[-1,10],xlabel='t[s]')
    ax2.legend(loc='upper right')
    
    plt.show()
    plt.savefig(os.path.join(img_path,sta,sta+'_Hkmap.png'))
    
    with open(os.path.join(img_path,sta,sta+'_Hk_TF.txt'),'w+') as f:
        f.write('[Hkpara]\n')
        f.write('station = %s\n'%sta)
        f.write('Vp = %.2f\n'%Vp)
        f.write('samprate = %.2f\n'%samprate)
        f.write('Hlimit[km] = %.1f-%.1f\n'%(Hlim[0],Hlim[1]))
        f.write('klimit = %.2f-%.2f\n'%(klim[0],klim[1]))
        f.write('H_search_range[km] = %.1f-%.1f\n'%(Hsearch[0],Hsearch[1]))
        f.write('k_search_range = %.2f-%.2f\n'%(ksearch[0],ksearch[1]))
        f.write('TFtype = %d\n'%TFtype)
        f.write('freq = %.2f,%.2f\n'%(freq[0],freq[1]))
        f.write('[Result]\n')
        f.write('H = %.2f\n'%Hfinal) 
        f.write('k = %.2f\n'%kfinal)
        if len(sigk) != 0 and len(sigh) != 0:
            f.write('Hstd = %.2f\n'%sigh[0])
            f.write('kstd = %.2f'%sigk[0])
        f.close()

def Hk_singletrace(data,rayp,Vp=6.5,samprate=10,Hlim=[4,8],Vslim=[3.,4.2],\
                   Hsearch=[4.5,7.5],Vssearch=[3,4.2],TFtype=1,freq=[0.1,2],RFwin=[10,30]):
    dH, dVs = 0.05, 0.05
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
    RF_time = np.linspace(-RFwin[0],RFwin[1],num=len(data),endpoint=True)

    Pamp = np.abs(data).max()
    RF_nomoho = propmat.rf_moho(0,0,rayp)
    for Vsid in range(len(Vsbg)):
        Vs = Vsbg[Vsid]
        for hid in range(len(Hbg)):
            h = Hbg[hid]
            print(h,Vs)
            try:
                RF_moho = propmat.rf_moho(h,Vs,rayp)
                # normalization
                tf = propmat.tf_cal(RF_nomoho/np.abs(RF_nomoho).max()*Pamp,RF_moho/np.abs(RF_moho).max()*Pamp)
                #tf = propmat.tf_cal(RF_nomoho,RF_moho)
                rf_clean = ifft(fft(data)*tf).real
                rf_clean_seg = rf_clean[np.where((RF_time>1.5)&(RF_time<4))]
                s[Vsid,hid] += np.linalg.norm(rf_clean_seg)
            except:
                print('synthetic failed')
                s[Vsid,hid] += np.max(s)
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
    
    # clean multiples using best (H,Vs) pair
    RF_moho = propmat.rf_moho(Hfinal,Vsfinal,rayp)
    tf = propmat.tf_cal(RF_nomoho/np.abs(RF_nomoho).max()*Pamp,RF_moho/np.abs(RF_moho).max()*Pamp)
    #tf = propmat.tf_cal(RF_nomoho,RF_moho)
    rf_clean = ifft(fft(data)*tf).real
    
    print('H=%.2f[km], k=%.2f'%(Hfinal,Vsfinal))

    # uncertainty
    sigs = np.std(s.reshape(1,-1)[0])
    pdVs2 = np.diff(s[i-1:i+2,j],2)/np.power(dVs,2)
    pdh2 = np.diff(s[i,j-1:j+2],2)/np.power(dH,2)
    sigVs = np.sqrt(2*sigs/np.abs(pdVs2))
    sigh = np.sqrt(2*sigs/np.abs(pdh2))
    
    fig = plt.figure(figsize=(6,8))
    ax1 = fig.add_axes([0.05,0.05,0.85,0.6])
    ax2 = fig.add_axes([0.05,0.7,0.85,0.25])
    ax3 = fig.add_axes([0.92,0.05,0.02,0.6]) # colorbar
    
    cb = ax1.pcolor(Hbg,Vsbg,s,cmap=plt.cm.Blues)
    ax1.add_patch(patches.Rectangle((Hsearch[0],Vssearch[0]),Hsearch[-1]-Hsearch[0],\
        Vssearch[-1]-Vssearch[0],fill=False,linewidth=1,linestyle='dashed'))
    ax1.vlines(Hfinal,Vsbg[0],Vsbg[-1],color='r',linestyle='dashed')
    ax1.hlines(Vsfinal,Hbg[0],Hbg[-1],color='r',linestyle='dashed')
    try:
        ex, ey = elipab([Hfinal,Vsfinal],sigh,sigVs)
        ax1.plot(ex,ey,'Vs')
        print('Hstd=%.2f[km], Vsstd=%.2f'%(sigh[0],sigVs[0]))
    except:
        pass
    fig.colorbar(cb,cax=ax3,label='RF norm')
    ax1.set(title='',xlabel='H[km]',ylabel='Vs[km/s]',xlim=Hlim,ylim=Vslim)
    ax2.plot(RF_time,data,label='original')
    ax2.plot(RF_time,rf_clean,label='cleaned')
    ax2.set(xlim=[-1,10],xlabel='t[s]')
    ax2.legend(loc='upper right')
    plt.show()
    plt.savefig('test_Hk.png')

def read_input(patch):
    sta = patch.entry01.get()
    Vp = float(patch.entry03.get())
    samprate = float(patch.entry05.get())
    Hlim = [float(patch.entry11.get()),float(patch.entry13.get())]
    klim = [float(patch.entry21.get()),float(patch.entry23.get())]
    Hsearch = [float(patch.entry31.get()),float(patch.entry33.get())]
    ksearch = [float(patch.entry41.get()),float(patch.entry43.get())]
    TFtype = int(patch.entry51.get())
    freq = [float(patch.entry53.get()),float(patch.entry55.get())]
    head = patch.entry61.get()
    config.read(head)
    prefix = head.replace(os.path.basename(head),'')
    out_path = os.path.join(prefix,config.get('path', 'out_path'))
    img_path = os.path.join(out_path,config.get('path', 'image_dir'))
    RF_path = os.path.join(out_path,config.get('path', 'RF_dir'))
    Hkstack(sta,Vp,samprate,Hlim,klim,Hsearch,ksearch,TFtype,freq,RF_path,img_path)

def loadHkparafile(patch):
    sta = patch.entry01.get()
    head = patch.entry61.get()
    config.read(head)
    out_path = config.get('path', 'out_path')
    img_path = os.path.join(out_path,config.get('path', 'image_dir'))
    Hkparafile = os.path.join(img_path,sta,sta+'_Hk_TF.txt')
    if not os.path.isfile(Hkparafile):
        print("Hkparafile doesn't exist!")
    try:
        config.read(Hkparafile)
        Vp = config.get('Hkpara','Vp')
        samprate = config.get('Hkpara','samprate')
        Hlim = config.get('Hkpara','Hlimit[km]').split('-')
        klim = config.get('Hkpara','klimit').split('-')
        Hsearch = config.get('Hkpara','H_search_range[km]').split('-')
        ksearch = config.get('Hkpara','k_search_range').split('-')
        TFtype = config.get('Hkpara','TFtype')
        freq = config.get('Hkpara','freq').split(',')
        
        patch.entry03.delete(0,'end'); patch.entry03.insert(0,Vp)
        patch.entry05.delete(0,'end'); patch.entry05.insert(0,samprate)
        patch.entry11.delete(0,'end'); patch.entry11.insert(0,Hlim[0])
        patch.entry13.delete(0,'end'); patch.entry13.insert(0,Hlim[1])
        patch.entry21.delete(0,'end'); patch.entry21.insert(0,klim[0])
        patch.entry23.delete(0,'end'); patch.entry23.insert(0,klim[1])
        patch.entry31.delete(0,'end'); patch.entry31.insert(0,Hsearch[0])
        patch.entry33.delete(0,'end'); patch.entry33.insert(0,Hsearch[1])
        patch.entry41.delete(0,'end'); patch.entry41.insert(0,ksearch[0])
        patch.entry43.delete(0,'end'); patch.entry43.insert(0,ksearch[1])
        patch.entry51.delete(0,'end'); patch.entry51.insert(0,TFtype)
        patch.entry53.delete(0,'end'); patch.entry53.insert(0,freq[0])
        patch.entry55.delete(0,'end'); patch.entry55.insert(0,freq[1])
    except:
        print('Load Hkparafile error!')

def HkGUI():
    win = tkinter.Tk()
    win.title('H-k stacking params')
    class Patch(object):
        pass
    patch = Patch()
    patch.label00 = tkinter.Label(win,text='sta:')
    patch.label02 = tkinter.Label(win,text='Vp:')
    patch.label04 = tkinter.Label(win,text='km/s   fs:')
    patch.label10 = tkinter.Label(win,text='H limit:')
    patch.label12 = tkinter.Label(win,text='to')
    patch.label14 = tkinter.Label(win,text='km')
    patch.label20 = tkinter.Label(win,text='k limit:')
    patch.label22 = tkinter.Label(win,text='to')
    patch.label30 = tkinter.Label(win,text='H search range:')
    patch.label32 = tkinter.Label(win,text='to')
    patch.label34 = tkinter.Label(win,text='km')
    patch.label40 = tkinter.Label(win,text='k search range:')
    patch.label42 = tkinter.Label(win,text='to')
    # 1-only Moho 2-only crustal multiples PpPs, 2P2S
    patch.label50 = tkinter.Label(win,text='TF type(1/2):')
    patch.label52 = tkinter.Label(win,text='Freq(Hz) lo:')
    patch.label54 = tkinter.Label(win,text='hi:')
    patch.label60 = tkinter.Label(win,text='RFparafile')
    patch.label62 = tkinter.Label(win,text='Hkparafile')
    patch.entry01 = tkinter.Entry(win,width=10)
    patch.entry03 = tkinter.Entry(win,width=10)
    patch.entry05 = tkinter.Entry(win,width=10)
    patch.entry11 = tkinter.Entry(win,width=10)
    patch.entry13 = tkinter.Entry(win,width=10)
    patch.entry21 = tkinter.Entry(win,width=10)
    patch.entry23 = tkinter.Entry(win,width=10)
    patch.entry31 = tkinter.Entry(win,width=10)
    patch.entry33 = tkinter.Entry(win,width=10)
    patch.entry41 = tkinter.Entry(win,width=10)
    patch.entry43 = tkinter.Entry(win,width=10)
    patch.entry51 = tkinter.Entry(win,width=10)
    patch.entry53 = tkinter.Entry(win,width=10)
    patch.entry55 = tkinter.Entry(win,width=10)
    patch.entry61 = tkinter.Entry(win,width=20)

    patch.label00.grid(row=0,column=0)
    patch.label02.grid(row=0,column=2)
    patch.label04.grid(row=0,column=4)
    patch.label10.grid(row=1,column=0)
    patch.label12.grid(row=1,column=2)
    patch.label14.grid(row=1,column=4)
    patch.label20.grid(row=2,column=0)
    patch.label22.grid(row=2,column=2)
    patch.label30.grid(row=3,column=0)
    patch.label32.grid(row=3,column=2)
    patch.label34.grid(row=3,column=4)
    patch.label40.grid(row=4,column=0)
    patch.label42.grid(row=4,column=2)
    patch.label50.grid(row=5,column=0)
    patch.label52.grid(row=5,column=2)
    patch.label54.grid(row=5,column=4)
    patch.label60.grid(row=6,column=0)
    patch.entry01.grid(row=0,column=1)
    patch.entry03.grid(row=0,column=3)
    patch.entry05.grid(row=0,column=5)
    patch.entry11.grid(row=1,column=1)
    patch.entry13.grid(row=1,column=3)
    patch.entry21.grid(row=2,column=1)
    patch.entry23.grid(row=2,column=3)
    patch.entry31.grid(row=3,column=1)
    patch.entry33.grid(row=3,column=3)
    patch.entry41.grid(row=4,column=1)
    patch.entry43.grid(row=4,column=3)
    patch.entry51.grid(row=5,column=1)
    patch.entry53.grid(row=5,column=3)
    patch.entry55.grid(row=5,column=5)
    patch.entry61.grid(row=6,column=1)

    patch.entry03.insert(0,'6.5')
    patch.entry05.insert(0,'10')
    patch.entry11.insert(0,'4')
    patch.entry13.insert(0,'8')
    patch.entry21.insert(0,'1.5')
    patch.entry23.insert(0,'2')
    patch.entry31.insert(0,'4.5')
    patch.entry33.insert(0,'7.5')
    patch.entry41.insert(0,'1.5')
    patch.entry43.insert(0,'2')
    patch.entry51.insert(0,'1')
    patch.entry53.insert(0,'0.1')
    patch.entry55.insert(0,'2')
    patch.entry61.insert(0,'../../paraRF.cfg')

    button63 = tkinter.Button(win,text="load Hkparafile",command=lambda: loadHkparafile(patch))
    button7 = tkinter.Button(win,text="search begin",command=lambda: read_input(patch))
    button63.grid(row=6,column=3)
    button7.grid(row=7)
    win.mainloop()

if __name__ == '__main__':
    obs_file = 'PCAstack/RF50_0.10to2.00/CC04_PCA.dat'
    data = np.loadtxt(obs_file)[:,1][::5]
    Hk_singletrace(data,0.06)