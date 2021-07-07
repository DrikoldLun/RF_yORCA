import os, configparser
#stalst = os.listdir('../RF_data/event_data_Pcorrected0.1')
parafile = 'paraRF.cfg'
config = configparser.ConfigParser()
config.read(parafile)
event_path = config.get('path','event_path')
#stalst = os.listdir(event_path)
stalst1 = ['CC06','EE04','CC11','WC05','WW03'] #0.8hz noise
stalst2 = ['CC02','CC04','CC05','CC07','CC08','WC03','WC04','WC02'] #normal stations
stalst = stalst1 + stalst2
#stalst = ['EE02'] #0.8hz noise
#stalst1 = ['EC05','EE01','WW01','WW04','CC03','EC01','CC12','CC09'] #0.8hz noise
#stalst2 = ['CC01','WC02'] #normal
#stalst = stalst1+stalst2

if not os.path.exists('log'):
    os.mkdir('log')
for sta in stalst:
    #os.system('rm -rf ../RF_data/RF50_0.10to2.00/%s'%sta)
    runlog = 'log/run_'+sta+'.log'
    runerr = 'log/run_'+sta+'.err'
    os.system('nohup python -u RFcal.py -S'+sta+' -O paraRF.cfg 1> '+runlog+' 2> '+runerr+' &')
    #os.system('nohup python -u rtzcal_event.py -S'+sta+' -O paraRF.cfg 1> '+runlog+' 2> '+runerr+' &')
