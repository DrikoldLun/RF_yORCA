import numpy as np

def rms(data):
    return np.sqrt(np.mean(np.power(data,2)))

def snr(data,fs,t_boundary,time_before,time_after):
    tb = t_boundary-time_before
    te = t_boundary+time_after
    nb = int(tb*fs-1)
    nc = int(t_boundary*fs-1)
    ne = int(te*fs-1)

    n_win = data[nb:nc]
    s_win = data[nc+1:ne+1]

    return max(np.abs(s_win))/rms(n_win)
