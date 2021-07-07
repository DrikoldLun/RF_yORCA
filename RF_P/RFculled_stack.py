import os, glob
import numpy as np
import matplotlib.pylab as plt

def plotRF(ax, data, rf_time, offsety, norm=1):
        y1, y2 = data.copy()*norm, data.copy()*norm
        thre = 0.03*np.max(np.abs(data*norm))
        y1[np.where(data<=thre)[0]] = thre
        y2[np.where(data>-thre)[0]] = -thre
        ax.fill_between(rf_time, thre+offsety, y1+offsety, facecolor='red',alpha=0.5)
        ax.fill_between(rf_time, -thre+offsety, y2+offsety, facecolor='blue',alpha=0.5)
        ax.plot(rf_time,data*norm+offsety,'k-',linewidth=0.5,alpha=0.5)

dat_path = '../RF_data/CCP/RF50_0.10to0.50/RFdat_culled'
#dat_rmul_path = '../RF_data/CCP/RF50_0.10to2.00/RFdat_rmul'
stalst = os.listdir(dat_path)
stalst.sort()
RFtime = np.linspace(-10,30,num=401,endpoint=True)
stackRF = np.zeros([len(stalst),len(RFtime)])
stackall = np.zeros(len(RFtime))
pamp = []

for i in range(len(stalst)):
    sta_RF_lst = glob.glob(os.path.join(dat_path,stalst[i],'*_PRF.R'))
    for rf_file in sta_RF_lst:
        stackRF[i] += np.loadtxt(rf_file)
    stackRF[i] /= len(sta_RF_lst)
    pamp.append(np.abs(stackRF[i]).max())
    stackall += stackRF[i]
stackall /= len(stalst)

norm = 0.8/np.median(pamp)

fig = plt.figure(figsize=(4,6))
ax = fig.add_subplot(111)

for i in range(len(stalst)):
    plotRF(ax,stackRF[i],RFtime,i,norm)
plotRF(ax,stackall,RFtime,len(stalst)+1,norm)
ax.set(xlabel='t[s]',xlim=[-1,10],ylabel='sta',yticks=np.arange(0,len(stalst)).tolist()+[len(stalst)+1],yticklabels=stalst+['stack'])
ax.grid()
plt.savefig('RF_culled_stack2hz.png')
plt.show()

'''
dat_path = '../RF_data/CCP/RF50_0.10to2.00/RFdat_culled'
dat_rmul_path = '../RF_data/CCP/RF50_0.10to2.00/RFdat_rmul'
stalst = os.listdir(dat_path)
stalst.sort()
RFtime = np.linspace(-10,30,num=401,endpoint=True)
stackRF = np.zeros([len(stalst),len(RFtime)])
stackall = np.zeros(len(RFtime))
stackRF_rmul = np.zeros([len(stalst),len(RFtime)])
stackall_rmul = np.zeros(len(RFtime))
pamp = [] 

for i in range(len(stalst)):
    sta_RF_lst = glob.glob(os.path.join(dat_path,stalst[i],'*_PRF.R'))
    for rf_file in sta_RF_lst:
        stackRF[i] += np.loadtxt(rf_file)
    stackRF[i] /= len(sta_RF_lst)
    pamp.append(np.abs(stackRF[i]).max())
    stackall += stackRF[i]
stackall /= len(stalst)

for i in range(len(stalst)):
    sta_RF_rmul_lst = glob.glob(os.path.join(dat_rmul_path,stalst[i],'*_PRF.R'))
    for rf_file in sta_RF_rmul_lst:
        stackRF_rmul[i] += np.loadtxt(rf_file)
    stackRF_rmul[i] /= len(sta_RF_rmul_lst)
    stackall_rmul += stackRF_rmul[i]
stackall_rmul /= len(stalst)

norm = 0.8/np.median(pamp)
stackRF *= norm
stackall *= norm
stackRF_rmul *= norm
stackall_rmul *= norm

fig = plt.figure(figsize=(4,6))
ax = fig.add_subplot(111)
for i in range(len(stalst)):
    ax.plot(RFtime,stackRF[i]+i,'b-')
    ax.plot(RFtime,stackRF_rmul[i]+i,'r-')

ax.plot(RFtime,stackall+len(stalst)+1,'b-',label='original')
ax.plot(RFtime,stackall_rmul+len(stalst)+1,'r-',label='crustal multiple removed')
ax.set(xlabel='t[s]',xlim=[-1,10],ylabel='sta',yticks=np.arange(0,len(stalst)).tolist()+[len(stalst)+1],yticklabels=stalst+['stack'])
ax.grid()
ax.legend(loc='upper right')
plt.savefig('rmul_stack.png')
plt.show()
'''
