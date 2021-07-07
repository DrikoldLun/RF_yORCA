import os
#stalst = os.listdir('../RF_data/event_data_Pcorrected0.1')
#parafile = 'paraRF.cfg'
#config = configparser.ConfigParser()
#config.read(parafile)
#event_path = config.get('path','event_path')
#stalst = os.listdir(event_path)
#stalst = ['CC06','EE04','CC07']
#stalst = ['CC04','CC05','WC03']
#stalst = ['CC05']
#stalst = ['CC02','CC04','CC05','CC07','CC08','WC03','WC04']
stalst = ['CC02','CC04','CC05','CC07','CC08','WC03','WC02']
#stalst = ['CC07','CC08']
if not os.path.exists('log'):
    os.mkdir('log')
for sta in stalst:
    runlog = 'log/run_'+sta+'_PCAstacktmp.log'
    runerr = 'log/run_'+sta+'_PCAstacktmp.err'
    os.system('nohup python -u propmat.py -S'+sta+' 1>> '+runlog+' 2>> '+runerr+' &')
    #os.system('nohup python -u checkrotate.py -S'+sta+' -L paraRF.cfg 1>> '+runlog+' 2>> '+runerr+' &')
