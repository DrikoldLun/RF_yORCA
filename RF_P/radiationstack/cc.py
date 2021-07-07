import numpy as np

def xcorr(x, y, lagmax=0, scale='none'):
    # Pad shorter array if signals are different lengths
    if x.size > y.size:
        pad_amount = x.size - y.size
        y = np.append(y, np.repeat(0, pad_amount))
    elif y.size > x.size:
        pad_amount = y.size - x.size
        x = np.append(x, np.repeat(0, pad_amount))

    corr = np.correlate(x, y, mode='full')  # scale = 'none'
    lags = np.arange(-(x.size - 1), x.size)

    if scale == 'biased':
        corr = corr / x.size
    elif scale == 'unbiased':
        corr /= (x.size - abs(lags))
    elif scale == 'coeff':
        corr /= np.sqrt(np.dot(x, x) * np.dot(y, y))

    if lagmax != 0:
        corr = corr[np.abs(lags)<=lagmax]
        lags = lags[np.abs(lags)<=lagmax]
    return corr, lags

def std_mccc(G,d):
    nsta = G.shape[1]
    nmat = G.shape[0]-1
    std = np.zeros(nsta)
    for i in range(nsta):
        dvcstd = np.zeros(nmat)
        for j in range(nmat):
            #if G[j,i] != 0: continue
            G_cut = np.vstack((G[:j],G[j+1:]))
            d_cut = np.hstack((d[:j],d[j+1:]))
            A = np.dot(G_cut.T,G_cut)
            Gg = np.dot(np.linalg.inv(A),G_cut.T)
            dcor = np.dot(Gg,d_cut)
            resid = d_cut[:-1]-np.dot(G_cut[:-1],dcor)
            dvcstd[j] = np.sqrt(np.sum(resid[G_cut[:-1,i]!=0]**2)/(np.sum(G_cut[:-1,i]!=0)-2))
        std[i] = np.std(dvcstd)
    return std

def mccc(data,dt,lagmax=1.,ifwt=False):
    # data = (npt x nsta) array of data
    # dt = 1/samplerate
    # 'lagmax' = max allowed lag (in seconds, default=1)
    # wt = 1D weighting array
    
    nsta = data.shape[1]
    nmat = nsta*(nsta-1)/2
    G = np.zeros([nmat,nsta])
    d = np.zeros(nmat)
    wvec = np.zeros(nmat)
    kk = 0

    for ii in range(nsta-1):
        for jj in range(ii+1,nsta):
            #[c,lags] = eng.xcorr(matlab.double(data[:,ii].tolist()),matlab.double(data[:,jj].tolist()),lagmax/dt,'unbiased',nargout=2)
            [c,lags] = xcorr(data[:,ii],data[:,jj],lagmax=int(lagmax/dt),scale='unbiased')
            pk = np.where(c==max(c))[0]
            G[kk,ii] = 1
            G[kk,jj] = -1
            d[kk] = lags[pk]*dt
            wvec[kk] = max(c)
            kk += 1
    
    G = np.vstack((G,np.ones(nsta)))
    d= np.hstack((d,np.zeros(1)))
            
    if not ifwt:
        wmat = np.identity(nmat+1)
    else:
        wmat = np.diag(np.hstack((wvec,np.ones(1))))
    GT = np.dot(G.T,wmat)
    A = np.dot(GT,G)
    Gg = np.dot(np.linalg.inv(A),GT)
    dcor = np.dot(Gg,d)
    
    resid = d[:kk]-np.dot(G[:kk,:],dcor)
    dcstd = np.sqrt(np.diag(np.linalg.inv(np.dot(G.T,G)))*np.var(resid))
    dvcstd = np.zeros(nsta)
    stack = np.zeros(len(data[:,0]))
    time_axis = np.arange(0,len(data[:,0]))*dt
    for ista in range(nsta):
        dvcstd[ista] = np.sqrt(np.sum(resid[G[:kk,ista]!=0]**2)/(np.sum(G[:kk,ista]!=0)-2))
        stack += np.interp(time_axis,time_axis-dcor[ista],data[:,ista])
    stack /= nsta
    ccmax = np.zeros(nsta)
    for ista in range(nsta):
        c, shift = xcorr(data[:,ista],stack,lagmax=int(lagmax/dt),scale='coeff')
        ccmax[ista] = max(c)

    rsstd = std_mccc(G,d)
    return dcor, dcstd, dvcstd, ccmax, rsstd

def ZP_mccc(staz, dataz, dtz, stap, datap, dtp, lagmax=1, ifwt=False):
    # stacomp = sta list for 1 comp
    # dcomp = cc between 2 traces for 1 comp (from mccc_matrice)
    sta = np.unique(staz+stap).tolist()
    nsta, nstaz, nstap =len(sta), len(staz), len(stap)
    nmatz, nmatp = nstaz*(nstaz-1)/2, nstap*(nstap-1)/2
    iz = [sta.index(tmp) for tmp in staz]
    ip = [sta.index(tmp) for tmp in stap]
    Gz = np.zeros([nmatz,nsta])
    Gp = np.zeros([nmatp,nsta])
    dz = np.zeros(nmatz)
    dp = np.zeros(nmatp)
    wtz = np.zeros(nmatz)
    wtp = np.zeros(nmatp)
    kk = 0
    for ii in range(nstaz-1):
        for jj in range(ii+1,nstaz):
            [c,lags] = xcorr(dataz[:,ii],dataz[:,jj],lagmax=int(lagmax/dtz),scale='unbiased')
            pk = np.where(c==max(c))[0][0]
            Gz[kk,iz[ii]] = 1
            Gz[kk,iz[jj]] = -1
            dz[kk] = lags[pk]*dtz
            wtz[kk] = max(c)
            kk += 1
    kk = 0
    for ii in range(nstap-1):
        for jj in range(ii+1,nstap):
            [c,lags] = xcorr(datap[:,ii],datap[:,jj],lagmax=int(lagmax/dtp),scale='unbiased')
            pk = np.where(c==max(c))[0][0]
            Gp[kk,ip[ii]] = 1
            Gp[kk,ip[jj]] = -1
            dp[kk] = lags[pk]*dtp
            wtp[kk] == max(c)
            kk += 1
    #print('Z')
    #print(wtz)
    #print('P')
    #print(wtp)
    G = np.vstack((Gz,Gp,np.ones(len(sta))))
    d = np.hstack((dz,dp,np.zeros(1)))
    
    wvec = np.hstack((wtz,wtp,np.ones(1)))
    
    if not ifwt:
        wmat = np.identity(nmatz+nmatp+1)
    else:
        wmat = np.diag(wvec)
        
    GT = np.dot(G.T,wmat)
    A = np.dot(GT,G)
    Gg = np.dot(np.linalg.inv(A),GT)
    dcor = np.dot(Gg,d)

    resid = d[:nmatz+nmatp]-np.dot(G[:nmatz+nmatp,:],dcor)
    dcstd = np.sqrt(np.diag(np.linalg.inv(np.dot(G.T,G)))*np.var(resid))
    dvcstd = np.zeros(nsta)
    for ista in range(nsta):
        dvcstd[ista] = np.sqrt(np.sum(resid[G[:nmatz+nmatp,ista]!=0]**2)/(np.sum(G[:nmatz+nmatp,ista]!=0)-2))
        
    rsstd = std_mccc(G,d)
    return sta, dcor, dcstd, dvcstd, rsstd
