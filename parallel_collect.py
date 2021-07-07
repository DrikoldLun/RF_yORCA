import os
#stalst = ['EE02']
stalst = ['CC01', 'CC09',  'CC12',  'EC05',  'EE01',  'WC02',  'WW01',  'WW04'] #valid drift
#stalst = ['CC03','EC01'] #invalid drift
#stalst = ['CC01',  'CC03', 'CC09',  'CC12',  'EC01',  'EC05',  'EE01',  'WC02',  'WW01',  'WW04'] #new12
if not os.path.exists('log'):
    os.mkdir('log')
mseed_path = '/home/lun/scratch-lun/OBSdata/yORCA_data/Mseed'
for sta in os.listdir(mseed_path):
    run_log = 'log/run_'+sta+'.log'
    run_err = 'log/run_'+sta+'.err'
    os.system('nohup python -u eventdata_collect.py -S'+sta+' 1> '+run_log+' 2> '+run_err+' &')
