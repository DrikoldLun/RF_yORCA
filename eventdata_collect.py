import sys, getopt
def Usage():
    print('python eventdata_collect.py -Sstation')
    sys.exit(1)
try:
    opts, args = getopt.getopt(sys.argv[1:],'S:')
except:
    Usage()
if opts == []:
    Usage()
for op, value in opts:
    if op == '-S':
        sta = value
    else:
        Usage()

import obspy,os
import numpy as np
from datetime import date
from obspy.taup import TauPyModel
#from obspy.clients.iris import Client
#import distaz

t_len = 3600
mseed_path = '../../../OBSdata/yORCA_data/Mseed'
culling_path = '../dataquality/PSD/hour_culling'
output_path = 'RF_data/event_data'
resp_path = '../data_file/RESP'

def preprocess(seis):
    seistmp = seis.copy()
    seistmp.detrend(type="demean")
    seistmp.detrend(type="linear")
    seistmp.taper(max_percentage=0.05)
    sta = seistmp.stats.station
    ch = seistmp.stats.channel
    respfile = os.path.join(resp_path,ch,'RESP.XX.'+sta+'..'+ch)
    pre_filt = (0.005, 0.01, 20.0, 25.0)
    seedresp = {'filename':respfile,
                'units':'VEL'}
    seistmp.simulate(pre_filt=pre_filt,seedresp=seedresp)
    return seistmp

with open('event_h.lst','r') as f:
    event_lst = f.readlines()
    f.close()

stlo, stla, stdp = {}, {}, {}
with open('../background/latlon_all.dat','r') as f:
    for line in f.readlines():
        line = line.replace('\n','').split()
        sta0, stlo0, stla0, stdp0 = line[0],float(line[1]), float(line[2]), float(line[3])
        stlo[sta0] = stlo0
        stla[sta0] = stla0
        stdp[sta0] = stdp0

'''
ev_num = 0
for e in range(len(event_lst)):
    event = event_lst[e].replace('\n','').split()
    hour_index = float(event[0])
    evlo = float(event[3])
    evla = float(event[4])
    evdp = float(event[5])
    tb = obspy.UTCDateTime(event[2])
    y, m, d = tb.year, tb.month, tb.day
    d_year = (date(y,m,d)-date(y,1,1)).days+1
    event_find = 0
    h_event = tb.hour
    min_event = tb.minute
    he_event = (tb+t_len).hour
    if he_event <= h_event:
        print('reject event '+str(e)+' for time discontinuity between 2 day records')
        continue
'''

#for sta in os.listdir(mseed_path):
def collect_sac(sta):
    culling_hour = []
    for i in range(3):
        CH = 'CH'+str(i)
        with open(os.path.join(culling_path,sta+'_'+CH+'_culling.dat'),'r') as f:
            for line in f.readlines():
                line = line.replace('\n','').split()
                culling_hour.append([int(tmp) for tmp in line[2:6]])
            f.close()
    evnum = 0
    for e in range(len(event_lst)):
        event = event_lst[e].replace('\n','').split()
        hour_index = float(event[0])
        evlo = float(event[3])
        evla = float(event[4])
        evdp = float(event[5])
        mag = float(event[1])

        tb = obspy.UTCDateTime(event[2])
        y, m, d = tb.year, tb.month, tb.day
        d_year = (date(y,m,d)-date(y,1,1)).days+1
        h_event = tb.hour
        min_event = tb.minute
        he_event = (tb+t_len).hour
        if he_event <= h_event:
            print('reject event '+str(e)+' for time discontinuity between 2 day records')
            continue

        tmp_path = os.path.join(mseed_path,sta,sta+'.*.'+str(y)+'.'+'%03d'%(d_year)+'*.msd')
        if len(os.popen('ls '+tmp_path).read()) == 0:
            continue
        '''
        client = Client()
        result = client.distaz(stla[sta],stlo[sta],evla,evlo)
        dis = result['distance']
        model = TauPyModel(model="iasp91")
        phase = ['P','S']
        arrivals = model.get_travel_times(source_depth_in_km=evdp,distance_in_degree=dis,phase_list=phase)
        if bodywave and len(arrivals) == 0:
            print('reject event '+str(e)+' for '+sta+'(no P,S arrival by TauP)')
            continue
        '''
        skip = False
        for h_inevent in range(h_event,he_event+1):
            if [y,m,d,h_inevent] in culling_hour:
                skip = True
                print('reject event '+str(e)+' for '+sta+'(event in culling_hour)')
                break
        if skip:
            continue
        new_ch = ['2','1','Z']
        for i in range(3):
            CH = 'CH'+str(i)
            msd = os.path.join(mseed_path,sta,'.'.join([sta,CH,str(y),'%03d'%d_year])+'*.msd')
            if len(os.popen('ls '+msd).read()) == 0:
                skip = True
                break
            st = obspy.read(msd)[0].copy()
            #st.data = st.data.astype('float64')
            fs = int(st.stats.sampling_rate)
            time_diff = int((tb-st.stats.starttime)*fs)
            ne = time_diff+fs*t_len
            if ne > len(st.data):
                ne = len(st.data)
            st.data = st.data[time_diff:ne]
            st.stats.starttime = tb
            try:
                st = preprocess(st)
            except:
                print("can't preprocess "+msd)
                skip = True
                break

            sta_path = os.path.join(output_path,sta)
            if not os.path.exists(sta_path):
                os.makedirs(sta_path)
            sac = os.path.join(sta_path,'%04d%02d%02d%02d%02d_%s.%1s'%(y,m,d,h_event,min_event,sta,new_ch[i]))
            st.write(sac,'SAC')
            try:
                st = obspy.read(sac)[0]
                st.stats.sac.evdp = evdp
                st.stats.sac.evla = evla
                st.stats.sac.evlo = evlo
                st.stats.sac.stla = stla[sta]
                st.stats.sac.stlo = stlo[sta]
                st.stats.sac.stdp = stdp[sta]
                st.stats.sac.mag = mag
                st.write(sac,'SAC')
            except:
                print("can't read "+sac)
                os.system('rm '+sac)
                skip = True
                break
            print(sac)

        if skip:
            continue
        evnum += 1
    print('%d events found in %s'%(evnum,sta))

stalst = os.listdir(mseed_path)
if sta in stalst:
    collect_sac(sta)
else:
    print("Station %s doesn't exit!"%sta)
