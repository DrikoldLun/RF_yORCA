#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 12:48:47 2019

@author: lun
"""

import sys, getopt, os, configparser, logging
# ------- Logger initiation ------- #
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
FORMAT = "line %(lineno)d - %(levelname)s: %(message)s"
formatter = logging.Formatter(FORMAT)
ch.setFormatter(formatter)
ch.setLevel(logging.INFO)
logger.addHandler(ch)

# ------- Read input ------- #
def Usage():
    msg = 'python rtzcal_event.py -Sstation [-O](reorientation?) [-F](frequency domain deconvolution) [-D](after sensorball drop) para.cfg'
    logger.info(msg)
    sys.exit(1)

try:
    opts, args = getopt.getopt(sys.argv[1:],'S:OFD')
except:
    Usage()
if opts == []:
    Usage()

reorient = False
fre_decov = False
afterdrop = False

for op, value in opts:
    if op == '-S':
        staname = value
    elif op == '-O':
        reorient = True
    elif op == '-F':
        fre_decov = True
    elif op == '-D':
        afterdrop = True
    else:
        Usage()

config = configparser.ConfigParser()
parafile = sys.argv[-1]
if os.path.isfile(parafile):
    config.read(parafile)

############################
# ------- Get Para ------- #
############################

event_path = config.get('path','event_path')
out_path = config.get('path','out_path')
rotate_path = config.get('path','rotate_dir')+'_event'
rotate_path = os.path.join(out_path,rotate_path)
RF_path = config.get('path','RF_dir')
RF_path = os.path.join(out_path,RF_path,staname)
image_path = config.get('path','image_dir')
image_path = os.path.join(out_path,image_path,staname)
for tmp_path in [rotate_path,RF_path,image_path]:
    if not os.path.exists(tmp_path): os.makedirs(tmp_path)
time_beforeP = config.getfloat('para','time_beforeP') 
time_afterP = config.getfloat('para','time_afterP')
gauss = config.getfloat('para','gauss')
freqmin = config.getfloat('para','freqmin')
freqmax = config.getfloat('para','freqmax')
dismin = config.getfloat('para','dismin')
dismax = config.getfloat('para','dismax')
magmin = config.getfloat('para','magmin')
samprate = config.getfloat('para','samprate')
noisegate = config.getfloat('para','noisegate')
RF_calcwinl = config.getfloat('para','RF_calcwinl')
RF_calcwinr = config.getfloat('para','RF_calcwinr')

X = 0
if reorient: 
    orient_file = config.get('path','orient_file')
    reorient_succeed = False
    with open(orient_file,'r') as f:
        for line in f.readlines():
            line = line.replace('\n','').split()
            if line[0] == staname:
                X = float(line[1])
                reorient_succeed = True
                break
        if not reorient_succeed:
            msg = "Specified orientation info of %s not found, original orientation used"%staname
            logger.warning(msg)       
        f.close()


if afterdrop:
    from obspy.core import UTCDateTime
    all_begin = UTCDateTime(2018,4,1)
    all_end = UTCDateTime(2019,6,1)
    drop_file = config.get('path','drop_file')
    drop_date = {}
    with open(drop_file,'r') as f:
        for line in f.readlines():
            line = line.replace('\n','').split()
            if line[1] == 'all':
                drop_date[line[0]] = all_end
            elif line[1] == 'no':
                drop_date[line[0]] = all_begin
            else:
                y,m,d = [int(tmp) for tmp in line[1].split('/')]
                drop_date[line[0]] = UTCDateTime(y,m,d)

###########################
# ----- Assign file ----- #
###########################

sys.path.append('RFscript')
from RFscript import decov, geo
import glob, obspy
import numpy as np
from obspy.taup import TauPyModel

#from obspy.clients.iris import Client
import distaz
from obspy.signal.rotate import rotate_ne_rt
from obspy.io.sac import SACTrace
import matplotlib.pylab as plt
model = TauPyModel(model="iasp91")

Z = glob.glob(os.path.join(event_path,staname,'*.Z'))
for z in Z:
    prefix = z[:-2].split('/')[-1]
    msg = "Start processing %s"%prefix
    logger.info(msg)
    try:
        st_z = obspy.read(z)[0].copy()
        st_1 = obspy.read(os.path.join(event_path,staname,prefix+'.1'))[0].copy()
        st_2 = obspy.read(os.path.join(event_path,staname,prefix+'.2'))[0].copy()
    except:
        msg = "channels of are not complete, continue"
        logger.info(msg)
        continue

###########################
# ----- Process RFs ----- #
###########################

# ------- Extract info ------- #
    tb = st_z.stats.starttime

    if afterdrop and tb < drop_date[staname]:
        msg = "evttime %s is eailier than sensorball drop date %s, continue"%(tb,drop_date[staname])
        logger.info(msg)
        continue
    
    stats = st_z.stats.sac
    result = distaz.DistAz(stats.stla,stats.stlo,stats.evla,stats.evlo)
    '''
    client = Client()
    result = client.distaz(stats.stla,stats.stlo,stats.evla,stats.evlo)
    dist = result['distancemeters']/1000.
    gcarc = result['distance']
    az = result['azimuth']
    baz = result['backazimuth']
    '''
    gcarc = result.getDelta()
    az = result.getAz()
    baz = result.getBaz()
    dist = result.degreesToKilometers(gcarc)

    if stats.mag < magmin:
        msg = "magnitude %.1f is small than threshold value %.1f, continue"%(stats.mag,magmin)
        logger.info(msg)
        continue

    if not dismin < gcarc < dismax:
        msg = "distance %.1f is not within range (%.1f,%.1f), continue"%(gcarc,dismin,dismax)
        logger.info(msg)
        continue

    try:
        arrivals = model.get_travel_times(source_depth_in_km=stats.evdp,distance_in_degree=gcarc,phase_list=['P'])
        t_P = arrivals[0].time
        rayp = geo.srad2skm(arrivals[0].ray_param)

    except:
        msg = "no theoretical P arrival, continue"
        logger.info(msg)
        continue

# ------- Rotate 12Z to RTZ ------- #
    
    if not len(st_1) == len(st_2) == len(st_z):
        msg = "Unequal record length among components"
        logger.error(msg)
        sys.exit(1)
    baz_used = np.mod(baz-X,360)
    data_r, data_t = rotate_ne_rt(st_1.data,st_2.data,baz_used)
    st_1.data, st_2.data = data_r, data_t
    st_r, st_t = st_1, st_2

# ------- head value assignment ------- #

    snr = {}
    for st_tmp,suffix in zip([st_r,st_t,st_z],['.R','.T','.Z']):
        st_tmp.stats.channel = suffix[1]
        st_tmp.stats.sac.dist = dist
        st_tmp.stats.sac.gcarc = gcarc
        st_tmp.stats.sac.az = az
        st_tmp.stats.sac.baz = baz
        st_tmp.stats.sac.user1 = freqmin
        st_tmp.stats.sac.user2 = freqmax
        st_tmp.stats.sac.user3 = RF_calcwinl
        st_tmp.stats.sac.user4 = rayp
        st_tmp.detrend(type='constant')
        st_tmp.detrend(type='linear')
        '''
        st_tmp.filter('bandpass',freqmin=freqmin,freqmax=freqmax,zerophase=True)
        fs = st_tmp.stats.sampling_rate
        if np.mod(fs,samprate) != 0:
            msg = "Designed samprate %d is not the factor of orignal samprate %d for %s"%(samprate,fs,prefix+suffix)
            logger.warning(msg)
        st_tmp.decimate(factor=int(fs//samprate),no_filter=True,strict_length=True)
        '''
# ------- Filtering and culling based on snr ------- #
        evt_path = os.path.join(rotate_path,prefix.split('_')[0])
        if not os.path.exists(evt_path):
            os.makedirs(evt_path)
        #st_tmp.write(os.path.join(rotate_path,prefix+suffix),'SAC')
        #st_tmp.filter('bandpass',freqmin=freqmin,freqmax=freqmax,zerophase=True)
        fs = st_tmp.stats.sampling_rate
        if np.mod(fs,samprate) != 0:
            msg = "Designed samprate %d is not the factor of orignal samprate %d for %s"%(samprate,fs,prefix+suffix)
            logger.warning(msg)
        st_tmp.decimate(factor=int(fs//samprate),no_filter=True,strict_length=True)
        #snr[suffix[1]] = geo.snr(st_tmp.data, samprate, t_S-2, 50, 15)
        #st_tmp.stats.sac.user5 = snr[suffix[1]]

    #if not [snr['R'],snr['Z']] > [noisegate]*2:
    #    msg = "snr not high enough, continue"
    #    logger.info(msg)
    #    continue

# ------- windowing ------- #
    
        winb, wine = int((t_P - RF_calcwinl)*samprate), int((t_P + RF_calcwinr)*samprate)
        st_tmp.data = st_tmp.data[winb:wine]
        st_tmp.stats.starttime += t_P - RF_calcwinl
        stsac = SACTrace.from_obspy_trace(st_tmp)
        stsac.reftime = st_tmp.stats.starttime + RF_calcwinl
        st_tmp = stsac.to_obspy_trace()
        st_tmp.write(os.path.join(evt_path,prefix+suffix),'SAC')
    print(prefix+' finished')
