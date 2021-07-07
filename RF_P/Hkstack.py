import numpy as np
import configparser, os, glob, obspy, tkinter
import matplotlib.pylab as plt
import matplotlib.patches as patches

config = configparser.RawConfigParser()

def elipab(o,lx,ly):
    theta = np.linspace(0,2*np.pi,20,endpoint=True)
    x = lx*np.cos(theta)+o[0]
    y = ly*np.sin(theta)+o[1]
    return x,y

def Hkstack(sta,Vp,Hlim,klim,Hsearch,ksearch,weight,RF_path,img_path):
    dH, dk = 0.1, 0.1
    nH = (Hlim[1]-Hlim[0])/dH+1
    nk = (klim[1]-klim[0])/dk+1
    Hbg = np.linspace(Hlim[0],Hlim[1],nH,endpoint=True)
    kbg = np.linspace(klim[0],klim[1],nk,endpoint=True)
    searchHb = int((Hsearch[0]-Hbg[0])/dH)
    searchHe = int((Hsearch[1]-Hbg[0])/dH)
    searchkb = int((ksearch[0]-kbg[0])/dk)+1
    searchke = int((ksearch[1]-kbg[0])/dk)+1
    s = np.zeros([len(kbg),len(Hbg)])
    s1, s2, s3 = s.copy(), s.copy(), s.copy()
    
    path = os.path.join(RF_path,sta)
    st = obspy.read(os.path.join(path,'*_PRF.R'))
    for kid in range(len(kbg)):
        k = kbg[kid]
        for hid in range(len(Hbg)):
            h = Hbg[hid]
            for rf in st:
                fs = rf.stats.sampling_rate
                tb = rf.stats.sac.b
                rayp = rf.stats.sac.user4
                Vs = Vp/k
                etas = np.sqrt(np.power(1./Vs,2)-np.power(rayp,2))
                etap = np.sqrt(np.power(1./Vp,2)-np.power(rayp,2))
                t1 = int((h*(etas-etap)-tb)*fs)
                t2 = int((h*(etas+etap)-tb)*fs)
                t3 = int((2*h*etas-tb)*fs)
                s[kid,hid] += weight[0]*rf.data[t1]+weight[1]*rf.data[t2]-weight[2]*rf.data[t3]
                s1[kid,hid] += weight[0]*rf.data[t1]
                s2[kid,hid] += weight[1]*rf.data[t2]
                s3[kid,hid] -= weight[2]*rf.data[t3]
    s /= len(st)
    s1 /= len(st)
    s2 /= len(st)
    s3 /= len(st)

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

    fig = plt.figure(figsize=(12,10))
    ax = []
    for i,stack,phase in zip(np.arange(0,4),[s,s1,s2,s3],['overall','Ps','PpPs','PsPs+PpSs']):
        ax.append(plt.axes([0.07+i%2*0.5,0.55-i//2*0.45,0.39,0.38]))
        cb = ax[i].pcolor(Hbg,kbg,stack,cmap=plt.cm.Blues)
        ax[i].add_patch(patches.Rectangle((Hsearch[0],ksearch[0]),Hsearch[-1]-Hsearch[0],ksearch[-1]-ksearch[0],fill=False,linewidth=1,linestyle='dashed'))
        ax[i].vlines(Hfinal,kbg[0],kbg[-1],color='r',linestyle='dashed')
        ax[i].hlines(kfinal,Hbg[0],Hbg[-1],color='r',linestyle='dashed')
        if i == 0:
            try:
                ex, ey = elipab([Hfinal,kfinal],sigh,sigk)
                ax[i].plot(ex,ey,'k')
                print('Hstd=%.2f[km], kstd=%.2f'%(sigh[0],sigk[0]))
            except:
                pass
        fig.colorbar(cb,ax=ax[i])
        ax[i].set(title=phase,xlabel='H[km]',ylabel='k',xlim=Hlim,ylim=klim)
    plt.savefig(os.path.join(img_path,sta,sta+'_Hkmap.png'))
    with open(os.path.join(img_path,sta,sta+'_Hk.txt'),'w+') as f:
        f.write('[Hkpara]\n')
        f.write('station = %s\n'%sta)
        f.write('Vp = %.2f\n'%Vp)
        f.write('Hlimit[km] = %.1f-%.1f\n'%(Hlim[0],Hlim[1]))
        f.write('klimit = %.2f-%.2f\n'%(klim[0],klim[1]))
        f.write('H_search_range[km] = %.1f-%.1f\n'%(Hsearch[0],Hsearch[1]))
        f.write('k_search_range = %.2f-%.2f\n'%(ksearch[0],ksearch[1]))
        f.write('Phase_weight = %.2f,%.2f,%.2f\n'%(weight[0],weight[1],weight[2]))
        f.write('[Result]\n')
        f.write('H = %.2f\n'%Hfinal) 
        f.write('k = %.2f\n'%kfinal)
        if len(sigk) != 0 and len(sigh) != 0:
            f.write('Hstd = %.2f\n'%sigh[0])
            f.write('kstd = %.2f'%sigk[0])
        f.close()
    plt.show()

def read_input():
    sta = entry01.get()
    Vp = float(entry03.get())
    Hlim = [float(entry11.get()),float(entry13.get())]
    klim = [float(entry21.get()),float(entry23.get())]
    Hsearch = [float(entry31.get()),float(entry33.get())]
    ksearch = [float(entry41.get()),float(entry43.get())]
    weight = [float(entry51.get()),float(entry53.get()),float(entry55.get())]
    head = entry61.get()
    config.read(head)
    out_path = config.get('path', 'out_path')
    img_path = os.path.join(out_path,config.get('path', 'image_dir'))
    RF_path = os.path.join(out_path,config.get('path', 'RF_dir'))
    Hkstack(sta,Vp,Hlim,klim,Hsearch,ksearch,weight,RF_path,img_path)

def loadHkparafile():
    sta = entry01.get()
    head = entry61.get()
    config.read(head)
    out_path = config.get('path', 'out_path')
    img_path = os.path.join(out_path,config.get('path', 'image_dir'))
    Hkparafile = os.path.join(img_path,sta,sta+'_Hk.txt')
    if not os.path.isfile(Hkparafile):
        print("Hkparafile doesn't exist!")
    try:
        config.read(Hkparafile)
        Vp = config.get('Hkpara','Vp')
        Hlim = config.get('Hkpara','Hlimit[km]').split('-')
        klim = config.get('Hkpara','klimit').split('-')
        Hsearch = config.get('Hkpara','H_search_range[km]').split('-')
        ksearch = config.get('Hkpara','k_search_range').split('-')
        weight = config.get('Hkpara','Phase_weight').split(',')
        
        entry03.delete(0,'end'); entry03.insert(0,Vp)
        entry11.delete(0,'end'); entry11.insert(0,Hlim[0])
        entry13.delete(0,'end'); entry13.insert(0,Hlim[1])
        entry21.delete(0,'end'); entry21.insert(0,klim[0])
        entry23.delete(0,'end'); entry23.insert(0,klim[1])
        entry31.delete(0,'end'); entry31.insert(0,Hsearch[0])
        entry33.delete(0,'end'); entry33.insert(0,Hsearch[1])
        entry41.delete(0,'end'); entry41.insert(0,ksearch[0])
        entry43.delete(0,'end'); entry43.insert(0,ksearch[1])
        entry51.delete(0,'end'); entry51.insert(0,weight[0])
        entry53.delete(0,'end'); entry53.insert(0,weight[1])
        entry55.delete(0,'end'); entry55.insert(0,weight[2])
    except:
        print('Load Hkparafile error!')

win=tkinter.Tk()
win.title('H-k stacking params')
label00 = tkinter.Label(win, text='sta:')
label02 = tkinter.Label(win, text='Vp:')
label04 = tkinter.Label(win, text='km/s')
label10 = tkinter.Label(win, text='H limit:')
label12 = tkinter.Label(win, text='to')
label14 = tkinter.Label(win, text='km')
label20 = tkinter.Label(win, text='k limit:')
label22 = tkinter.Label(win, text='to')
label30 = tkinter.Label(win, text='H search range:')
label32 = tkinter.Label(win, text='to')
label34 = tkinter.Label(win, text='km')
label40 = tkinter.Label(win, text='k search range:')
label42 = tkinter.Label(win, text='to')
label50 = tkinter.Label(win, text='Weight w1(Ps):')
label52 = tkinter.Label(win, text='w2(PpPs):')
label54 = tkinter.Label(win, text='w3(PpSs+PsPs):')
label60 = tkinter.Label(win, text='RFparafile')
label62 = tkinter.Label(win, text='Hkparafile')
entry01 = tkinter.Entry(win,width=10)
entry03 = tkinter.Entry(win,width=10)
entry11 = tkinter.Entry(win,width=10)
entry13 = tkinter.Entry(win,width=10)
entry21 = tkinter.Entry(win,width=10)
entry23 = tkinter.Entry(win,width=10)
entry31 = tkinter.Entry(win,width=10)
entry33 = tkinter.Entry(win,width=10)
entry41 = tkinter.Entry(win,width=10)
entry43 = tkinter.Entry(win,width=10)
entry51 = tkinter.Entry(win,width=10)
entry53 = tkinter.Entry(win,width=10)
entry55 = tkinter.Entry(win,width=10)
entry61 = tkinter.Entry(win,width=20)

label00.grid(row=0,column=0)
label02.grid(row=0,column=2)
label04.grid(row=0,column=4)
label10.grid(row=1,column=0)
label12.grid(row=1,column=2)
label14.grid(row=1,column=4)
label20.grid(row=2,column=0)
label22.grid(row=2,column=2)
label30.grid(row=3,column=0)
label32.grid(row=3,column=2)
label34.grid(row=3,column=4)
label40.grid(row=4,column=0)
label42.grid(row=4,column=2)
label50.grid(row=5,column=0)
label52.grid(row=5,column=2)
label54.grid(row=5,column=4)
label60.grid(row=6,column=0)
entry01.grid(row=0,column=1)
entry03.grid(row=0,column=3)
entry11.grid(row=1,column=1)
entry13.grid(row=1,column=3)
entry21.grid(row=2,column=1)
entry23.grid(row=2,column=3)
entry31.grid(row=3,column=1)
entry33.grid(row=3,column=3)
entry41.grid(row=4,column=1)
entry43.grid(row=4,column=3)
entry51.grid(row=5,column=1)
entry53.grid(row=5,column=3)
entry55.grid(row=5,column=5)
entry61.grid(row=6,column=1)

entry03.insert(0,'6')
entry11.insert(0,'10')
entry13.insert(0,'30')
entry21.insert(0,'1.5')
entry23.insert(0,'2')
entry31.insert(0,'10')
entry33.insert(0,'30')
entry41.insert(0,'1.5')
entry43.insert(0,'2')
entry51.insert(0,'0.7')
entry53.insert(0,'0.2')
entry55.insert(0,'0.1')
entry61.insert(0,'paraRF.cfg')

button63 = tkinter.Button(win,text="load Hkparafile",command=loadHkparafile)
button7 = tkinter.Button(win,text="search begin",command=read_input)
button63.grid(row=6,column=3)
button7.grid(row=7)
win.mainloop()
