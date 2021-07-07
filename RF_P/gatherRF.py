import configparser, os
import numpy as np
latlon_all = '../../background/latlon_all.dat'
stalst = np.loadtxt(latlon_all,usecols=0,dtype='str').tolist()
if __name__ == '__main__':
    config = configparser.RawConfigParser()
    parafile = 'paraRF.cfg'
    config.read(parafile)
    prefix = parafile.replace(os.path.basename(parafile),'')
    out_path = os.path.join(prefix,config.get('path', 'out_path'))
    img_path = os.path.join(out_path,config.get('path', 'image_dir'))
    RF_path = os.path.join(out_path,config.get('path', 'RF_dir'))
    evtlog_path = os.path.join(prefix,config.get('path', 'cull_evtlist'),config.get('path','RF_dir'))
    gather_path = os.path.join(RF_path,'all')
    if not os.path.exists(gather_path): os.mkdir(gather_path)
    for sta in stalst:
        goodrf_file = os.path.join(evtlog_path,'XX.%s_goodRF.lst'%sta)
        if not os.path.exists(goodrf_file):
            print('%s continue: not usable'%sta)
            continue
        sta_path = os.path.join(RF_path,sta)
        for evt_time in np.loadtxt(goodrf_file,dtype='str').tolist():
            os.system('cp %s/%s_%s_PRF.R %s/'%(sta_path,evt_time,sta,gather_path))
        print('%s finished'%sta)
