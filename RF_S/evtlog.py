import os, glob, obspy, configparser
import numpy as np

def evt_catalog(ch_path):
    evt_info = []
    # evttime, evlo, evla, evdp, mag, trnum
    evt_lst = os.listdir(ch_path)
    for evt in evt_lst:
        path_evt = os.path.join(ch_path,evt)
        sac_lst = glob.glob(os.path.join(path_evt,evt+'*'))
        if len(sac_lst) == 0: os.rmdir(path_evt)
        else:
            stats = obspy.read(sac_lst[0])[0].stats.sac
            evlo, evla, evdp, mag = stats.evlo, stats.evla, stats.evdp, stats.mag
            stanum = len(sac_lst)//3
            evt_info.append([evt,evlo,evla,evdp,mag,stanum])
    
    with open(os.path.join('evt.log'),'w+') as f:
        for evtline in evt_info:
            f.write(' '.join([str(tmp) for tmp in evtline])+'\n')
        f.close()

# read config
config = configparser.ConfigParser()
parafile = 'paraSRF.cfg'
if os.path.isfile(parafile): config.read(parafile)
else: Usage()
out_path = config.get('path','out_path')
ch_path = os.path.join(out_path,config.get('path','rotate_dir'))

evt_catalog(ch_path)
