import os, configparser
#stalst = os.listdir('../RF_data/event_data_Pcorrected0.1')
parafile = 'paraSRF.cfg'
config = configparser.ConfigParser()
config.read(parafile)
event_path = config.get('path','event_path')
stalst = os.listdir(event_path)
#stalst = ['CC06','EE04','CC07']
#stalst = ['CC02','CC04','CC05','CC07','CC08','WC03','WC04'] #normal stations
#stalst = ['EE02'] #0.8hz noise
#stalst1 = ['EC05','EE01','WW01','WW04','CC03','EC01','CC12'] #0.8hz noise
#stalst2 = ['CC01','WC02','CC09'] #normal
#stalst = stalst+stalst1+stalst2
if not os.path.exists('log'):
    os.mkdir('log')
for sta in stalst:
    runlog = 'log/run_'+sta+'.log'
    runerr = 'log/run_'+sta+'.err'
    os.system('nohup python -u RFcal.py -S'+sta+' -O paraSRF.cfg 1>> '+runlog+' 2>> '+runerr+' &')
    print(sta)
    #os.system('nohup python -u checkrotate.py -S'+sta+' -L paraRF.cfg 1>> '+runlog+' 2>> '+runerr+' &')
    #os.system('nohup python -u PCApick.py -S'+sta+' '+parafile+' 1>> '+runlog+' 2>> '+runerr+' &')
    #os.system('nohup python -u rtzcal_event.py -S'+sta+' -O paraSRF.cfg 1>> '+runlog+' 2>> '+runerr+' &')
