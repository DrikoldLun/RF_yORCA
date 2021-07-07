import os,sys,getopt,configparser,glob,obspy
import numpy as np

def Usage():
    print('python CCP_init.py -Cculling para.cfg')
    sys.exit(1)

try:
    opts, args = getopt.getopt(sys.argv[1:],'C')
except:
    Usage()

evtculling = False
for op, value in opts:
    if op == '-C':
        evtculling = True
    else:
        Usage()

if sys.argv[1:] == []: Usage()
parafile = sys.argv[-1]
if not os.path.isfile(parafile): Usage()

config = configparser.ConfigParser()
config.read(parafile)
out_path = config.get('path','out_path')
RF_path = os.path.join(out_path,config.get('path','RF_dir'))

if not evtculling:
    CCPdat_path = os.path.join(out_path,'CCP',os.path.basename(RF_path),'RFdat')

    if not os.path.exists(CCPdat_path): os.makedirs(CCPdat_path)
    for sta in os.listdir(RF_path):
        stadat_dir = os.path.join(CCPdat_path,sta)
        if not os.path.exists(stadat_dir): os.mkdir(stadat_dir)
        saclst = glob.glob(os.path.join(RF_path,sta,'*.R'))
        for sac in saclst:
            dat = os.path.join(stadat_dir,os.path.basename(sac))
            np.savetxt(dat,obspy.read(sac)[0].data,newline='\n')
        os.system('cp '+os.path.join(RF_path,sta,'*finallist.dat')+' '+stadat_dir)
else:
    evtlst_path = os.path.join(config.get('path', 'cull_evtlist'),config.get('path','RF_dir'))
    CCPdat_path = os.path.join(out_path,'CCP',config.get('path', 'RF_dir'),'RFdat_culled')
    if not os.path.exists(CCPdat_path): os.makedirs(CCPdat_path)
    for sta in os.listdir(RF_path):
        evtlist_file = os.path.join(evtlst_path,'XX.%s_goodRF.lst'%sta)
        if not os.path.exists(evtlist_file):
            print('%s continue'%sta)
            continue
        stadat_dir = os.path.join(CCPdat_path,sta)
        if not os.path.exists(stadat_dir): os.mkdir(stadat_dir)
        evtlist = np.loadtxt(evtlist_file).astype('int').astype('str')
        saclst = glob.glob(os.path.join(RF_path,sta,'*.R'))
        for sac in saclst:
            evt_time = os.path.basename(sac)[:12]
            if evt_time in evtlist:
                dat = os.path.join(stadat_dir,os.path.basename(sac))
                np.savetxt(dat,obspy.read(sac)[0].data[::5],newline='\n')
        finallist = os.path.join(evtlst_path,'XX.%sfinallist.dat'%sta)
        os.system('cp %s %s'%(finallist,stadat_dir))
