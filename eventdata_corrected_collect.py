import obspy, os, glob
from scipy.io import loadmat

raw_mat_path = '../dataquality/event/ATaCR/NOISETC_CI/DATA/datacache_prepro'
mat_path = '../dataquality/event/ATaCR/NOISETC_CI/DATA/NOISETC/CORRSEISAVTF0.1'
evt_file = 'event_h.lst'
sta_file = '../background/latlon.dat'
output_path = 'RF_data/event_data_Hcorrected0.1'

evt_info, sta_info = {}, {}
with open(evt_file,'r') as f:
    for line in f.readlines():
        line = line.replace('\n','').split()
        evlo, evla, evdp, mag = [float(tmp) for tmp in line[3:]+[line[1]]]
        tb = obspy.UTCDateTime(line[2])
        y, m, d, h, mi = tb.year, tb.month, tb.day, tb.hour, tb.minute
        evt_id = '%04d%02d%02d%02d%02d'%(y, m, d, h, mi)
        evt_info[evt_id] = [evlo, evla, evdp, mag, tb]
    f.close()

with open(sta_file,'r') as f:
    for line in f.readlines():
        line = line.replace('\n','').split()
        sta0, stlo0, stla0, stdp0 = line[0], float(line[1]), float(line[2]), float(line[3])
        sta_info[sta0] = [stlo0, stla0, stdp0]
    f.close()

for netsta in os.listdir(mat_path):
    sta_mat_path = os.path.join(mat_path,netsta)
    evnum = 0
    for mat_name in glob.glob(os.path.join(sta_mat_path,'*.mat')):
        mat = loadmat(mat_name)
        fs = 1./mat['corrected']['params'][0,0]['dt'][0,0][0,0]
        sta = mat['corrected']['params'][0,0]['station'][0,0][0]
        net = mat['corrected']['params'][0,0]['network'][0,0][0]
        evt_id = mat['corrected']['params'][0,0]['eventid'][0,0][0]
        correct_type = mat['corrseis']['label'][0]
        for i in range(len(correct_type)):
            if correct_type[i][0] == 'Z2-1': #'ZP-21':
                index = i
                break
        Z = mat['corrseis']['timeseries'][0,index].reshape(1,-1)[0]
        evlo, evla, evdp, mag, tb = evt_info[evt_id]
        stlo, stla, stdp = sta_info[sta]

        raw_mat = loadmat(os.path.join(raw_mat_path,evt_id,evt_id+'_'+net+'_'+sta+'.mat'))
        traces = raw_mat['traces'][0]
        for i in range(len(traces)):
            if traces[i]['channel'][0] == 'CH1':
                H1 = traces[i]['data'].reshape(1,-1)[0]
            elif traces[i]['channel'][0] == 'CH0':
                H2 = traces[i]['data'].reshape(1,-1)[0]
        for data, ch in zip([H1,H2,Z],['H1','H2','Z']):
            data = data[:int(3600*fs)]
            if len(data)/fs != 3600:
                print(evt_id+' data gap; continue')
                continue
            st = obspy.read()[0]
            st.stats.network = net
            st.stats.station = sta
            st.data = data
            st.stats.channel = ch
            st.stats.sampling_rate = fs
            st.stats.location = ''
            st.stats.starttime = tb
            sacdir = os.path.join(output_path,sta)
            if not os.path.exists(sacdir):
                os.makedirs(sacdir)
            sac = os.path.join(sacdir,evt_id+'_'+sta+'.'+ch[-1])
            st.write(sac,'SAC')
            try:
                st = obspy.read(sac)[0]
                st.stats.sac.evdp = evdp
                st.stats.sac.evla = evla
                st.stats.sac.evlo = evlo
                st.stats.sac.stla = stla
                st.stats.sac.stlo = stlo
                st.stats.sac.stdp = stdp
                st.stats.sac.mag = mag
                st.write(sac,'SAC')
            except:
                print("can't read "+sac)
                os.remove(sac)
                continue
        evnum += 1
    print('%d events found in %s'%(evnum,netsta))
        
