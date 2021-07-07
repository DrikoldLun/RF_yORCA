import sys, getopt, configparser, os
import numpy as np
import obspy
import plotrf
import matplotlib.pylab as plt

def Hdecomp(RFmat,baz,K=2):
    RFnum = RFmat.shape[0]
    G = np.zeros([RFnum,2*K+1]) #sinusoidal term
    for i in range(RFnum):
        G[i,0] = 1
        for j in range(1,K+1):
            G[i,2*j-1] = np.cos(j*np.deg2rad(baz[i]))
            G[i,2*j] =  np.sin(j*np.deg2rad(baz[i]))
    
    A = np.dot(G.T,G)
    Gg = np.dot(np.linalg.inv(A),G.T)
    H = np.dot(Gg,RFmat) #coefficient for sinusoidal terms
    return G, H

def Hsyn(G,H,order=[0]):
    col = []
    order = np.array(order)
    for i in order:
        if i == 0: 
            col.append(i)
        else:
            col.append(2*i-1)
            col.append(2*i)
    col = np.array(col).astype('int')
    G_syn = G[:,col]
    H_syn = H[col,:]
    RF_syn = np.dot(G_syn,H_syn)
    return RF_syn

def plotbaz(ax,RFmat,RFtime,baz):
    for i in range(len(baz)):
        plotrf.plotRF(ax,RFmat[i],RFtime,baz[i],norm=10)
    ax.set(xlim=[0,10],xticks=np.arange(0,10,2),xlabel='t[s]',ylim=[0,360],yticks=np.arange(0,360,30))
    ax.grid()

def Usage():
    print('python Hdecomp.py -Ssta -Cchannel parafile')
    sys.exit(1)
try:
    opts, args = getopt.getopt(sys.argv[1:],'S:C:')
except:
    Usage()
if opts == []:
    Usage()

for op, value in opts:
    if op == '-S':
        sta = value
    elif op == '-C':
        ch = value
    else:
        Usage()

head = sys.argv[-1]
config = configparser.RawConfigParser()
if os.path.isfile(head):
    config.read(head)
else:
    Usage()

out_path = config.get('path', 'out_path')
image_path = os.path.join(out_path,config.get('path', 'image_dir'),sta)
RF_path = os.path.join(out_path,config.get('path', 'RF_dir'),sta)
time_before = config.getint('para', 'time_beforeP') * -1
time_after = config.getint('para', 'time_afterP')
K = 2 #order

st = obspy.read(os.path.join(RF_path,'*.%s'%ch))
baz = [tr.stats.sac.baz for tr in st]
RFmat = np.zeros([len(st),st[0].stats.npts])
for i in range(len(st)):
    RFmat[i,:] = st[i].data[:]

G, H = Hdecomp(RFmat,baz,K=K)
RFtime = np.linspace(time_before,time_after,st[0].stats.npts)
fig = plt.figure(figsize=((K+4)*3,10))

ax1 = fig.add_subplot(1,K+4,1)
plotbaz(ax1,RFmat,RFtime,baz)
ax1.set_title('R')

ax2 = fig.add_subplot(1,K+4,2)
RF_H = Hsyn(G,H,order=np.arange(0,K+1))
plotbaz(ax2,RF_H,RFtime,baz)
ax2.set_title('H')

ax3 = fig.add_subplot(1,K+4,3)
plotbaz(ax3,RFmat-RF_H,RFtime,baz)
ax3.set_title('R-H')

for i in range(0,K+1):
    ax = fig.add_subplot(1,K+4,i+4)
    RF = Hsyn(G,H,order=[i])
    plotbaz(ax,RF,RFtime,baz)
    ax.set_title('A%d'%i)

np.savetxt('PCA/stack/RF50_0.10to2.00/Hconstant_%s.dat'%sta,np.vstack((RFtime,Hsyn(G,H,order=[0]))).T)
plt.suptitle(sta)
plt.savefig('Hdecomp/Hdecomp_%s.png'%sta)
plt.show()
