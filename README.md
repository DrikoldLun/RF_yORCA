# RF_yORCA

Please add RF_P/RFscript to python path

e.g. echo 'export PYTHONPATH=[your dir]/RF_yORCA/RF_P/RFscript:$PYTHONPATH' >> ~/.bashrc

## I. Event data collection (./)

script file: ***eventdata_collect.py***

Cut 1-hour seismic record since Earthquake origin time from daily mseed files

Usage: *python eventdata_collect.py -Sstation*

Input:  

(1) Line 25 Path for daily mseed files of seismic record

(2) Line 26 Path for culling records of stations

(3) Line 27 Output path (default: ***RF_data/event_data***)

(4) Line 28 Path for instrumental response files

(5) Line 44 Event list file  (default: ***event_h.lst***)

(6) Line 49 Latlon file for stations

Output: 1-hour sac files since event origin time for your station 

[your output address]/[station]/[eventtime]\_[station].[channel(1|2|Z)]

**Parallel mode:** specify **stalst** in ***parallel_collect.py***, then run '*python parallel_collect.py*'



## II. RF calculation (RF_P/)

(1) config file: ***paraRF.cfg***

**Important parameters:**

Path:

[event_path] The path for the event sac file collected in **I**

[out_path] The output path for rotated seismic data and RF data

[rotate_dir] the final path for rotated seismic data would be [out_path]/[rotate_dir]

[RF_dir] the final path for RF data would be [out_path]/[RF_dir] typically RF[samprate]_[freqmin]to[freqmax] like RF50_0.10to0.50

[orient_file] the station orientation file (default ***../result.txt***)

[drop_file] the hanging noise status file for stations

Parameters for RF calculation process:

[magmin] minimum magnitude threshold

[dismin] and [dismax] lower and upper distance bound

[noisegate] the minimum threshold for SNR

[gauss] the $\alpha$ parameter in deconvolution (I suggest 5 for investigating oceanic lithosphere)

[freqmin] and [freqmax] lower and upper frequency bound for filtering seismic record before deconvolution

[time_beforeP]  and [time_afterP] the final RF time would be [-time_beforeP, time_afterP]

[RF_calcwinlen] we would use the seismic record in [-RF_calcwinlen, RF_calcwinlen]s for deconvolution, the reference time is theoretical P arrival by TauP

[samprate] the sampling rate for RFs (I suggest 50 for ORCA case)



(2) RF calculation script: ***RFcal.py***

calculate RFs for a specific station

Usage: *python RFcal.py -Sstation [-O] (station orientation, needed for ORCA case) [-F] (frequency domain deconvolution) [-D] (only use records after sensorball drop) para.cfg*

For our case, typically run '*python RFcal.py -Sstation -O para.cfg*'

output: 

rotated seismic sac files in [out_path]/[rotate_dir]/[station]/[eventtime]_[station].[channel(R|T|Z)]

RF sac files in [out_path]/[RF_dir]/[station]/[eventtime]\_[station]\_PRF.[channel(R|T)]

**Parallel calculation mode**: specify **stalst** in ***RFcal_multi.py***, then run '*python RFcal_multi.py*'

Related python modules (deconvolution, distaz) are included in ***RFscript/***



## III. RF culling (RF_P/)

(1) config file: ***paraRF.cfg***

**Important parameters:**

[image_dir] the final image path would be [out_path]/[image_dir]

[cull_evtlist] path for the file recording culling result of RFs, please use the default 'PCA/evtlog'



(2) check the reliability of station orientation found

run 'python checkrotate.py -SCC06 paraRF.cfg'

output: figure in *PCA/checkrotate/CC06/*



(3) RF culling script: ***PCApickrf.py***

Usage:  python PCApickrf.py -Sstation -C[R/T]channel  para.cfg

typically use ' python PCApickrf.py -Sstation -CR  paraRF.cfg

As for the GUI program for RF culling, typically step:

ampcull->PCAcull->replot/deletebadRFs->manual cull(click the waveform)->Finish(save the culling result)

You can use 'Plot RFs' to plot the alignment of RFs in any stage

Output: 

In [cull_evtlist]/[RF_dir]/, XX.[station]\_badRF.lst, XX.[station]\_goodRF.lst, XX.[station]_PCA.dat(the stacked RF from retained RFs), XX.[station]finallist.dat (detailed information for retained RFs, please refer to Line 171 in ***PCApickrf.py*** for the meaning of different columns)

figures in [out_path]/[image_dir]/[station]/



(4) cull RFs for all stations concurrently

after finish culling for all individual stations

firstly run 'python gatherRF.py', output: all good RFs in [out_path]/[RF_dir]/all/

then run 'python PCApickrf.py -Sall -CR  paraRF.cfg' to cull these RFs concurrently



## IV. Hk stack (RF_P/)

run 'python Hkstack.py' and specify parameters in GUI, click 'load Hkparafile' to load existing Hkparafile, click 'search begin' to view Hk stack result

output: Hk stack image and parafile in [out_path]/[image_dir]/[station]/



## V. Harmonic decomposition

run 'python Hdecomp.py -Sstation -C[R|T]channel paraRF.cfg'

output: Harmonic decomposition figure in ***Hdecomp/***

H constant waveform  in PCA/stack/Hconstant_[station].dat



## VI. CCP migration

### 1. crustal multiple removal (RF_P/PROPMAT-master/matlab\_to\_propmat/)

you need to change Line 48 in 'writePROPMATexecfile_gauss.m' and Line 48 in 'writePROPMATexecfile.m' to specify your own working path

related script: ***propmat.py***

Usage: 'python propmat.py -Ssta'

Use Line 698 to get and save optimal crustal parameters (H and Vs) for crustal multiple cleanning

Output: HVs.dat (list for optimal crustal thickness H and Vs), HVs[sta].png (contour plot of the RF norm)

Use the code below to correct individual RFs

Output: RFs with crustal multiple removed in  ../../[out_path]/CCP/[RF_dir]/RFdat_rmul/[sta]/

Individual RF cleaning figure in  ../../[out_path]/[image_dir]/[station]/

Parallel mode: specify **stalst** in ***loop.py***, then run 'python loop.py'



### 2. water phase removal (RF_P/waterremove/)

related script: ***synth.py***

run *python synth.py -Ssta*

Input: RFs with crustal multiple removed in  ***../[out_path]/CCP/[RF_dir]/RFdat_rmul/[sta]/ (for 0.1-0.5Hz)***

 or ***../[out_path]/CCP/[RF_dir]/RFdat_culled/[sta]/ (for 0.1-2Hz)***

Output: 

RFs with crustal multiple and water phase removed (for 0.1-0.5Hz) or merely water phase removed (for 0.1-2Hz) in  ../[out_path]/CCP/[RF_dir]/RFdat_rwater/[sta]/

water phase removal effect figure RFwaterclean_[sta]\[freqmin]to[freqmax]\_culled.png



### 3. CCP migration (RF_P/CCP_0.5hz or RF_P/CCP_2hz)

in ***../***, run 

*'python CCP_init.py paraRF.cfg'* and *'python CCP_init.py -C paraRF.cfg'* respectively

output: all original RFs in ***../[out_path]/CCP/[RF_dir]/RFdat/*** and original RFs passing culling step in ***../[out_path]/CCP/[RF_dir]/RFdat_culled/***

There are 4 folders in  ***../[out_path]/CCP/[RF_dir]***: RFdat (all original RFs), RFdat_culled (original RFs passing culling step), RFdat_rmul (RFs culled with crustal multiple removed, not applicable for 0.1-2HZ), RFdat_rwater (RFs culled with water phase removed, and crustal multiple removed for 0.1-0.5Hz)

We have same operations in ***RF_P/CCP_0.5hz*** and ***RF_P/CCP_2hz***, they just differ in frequency band, CCP_0.5hz (0.1-0.5Hz), CCP_2hz(0.1-2Hz, crustal multiple removal doesn't apply due to timing)

Take ***RF_P/CCP_0.5hz*** for example:

Note: The name of all intermediate files must take the form 'filename\_[datatype(culled|rmul|rwater)]\_[stalst(newsta|oldsta)].xxx'

【1】firstly run 'PS_circlebin_makedata.m' to prepare ray tracing result, you need to specify

(1) Line 4, the station list

(2) Line 6, data path, you can use any of the four paths mentioned above

(3) Line 7, the velocity model, here 'IASP91.vel' is the crust-modified model, 'IASP91.velbak' is the original model

(4) Line 8, the coordinate of the profile

(5) Line 9, make sure 'shift' = 'time_beforeP'

(6) Line 115, the name of the RFdepth file, it would be loaded in the next script, like 'RFdepth_rwater_newsta'

Output: the RFdepth file containing ray tracing results

【2】then run 'Common_rectangularbin_stack_moho' to stack amplitude in each bins on the profile plane, you need to specify:

(1) Line 3, name of RFdepth file specified in the prior script

(2) Line 5, name of the output CCP stack datfile, this would be used in the plotting script

(3) Line 6, the coordinate of the profile

(4) Line 76, name of the file containing the projected station coordinate

output:

The .dat file containing stacked amplitude in each bins, like 'ccp_0.1-0.5hz_rwater_newsta.dat'

format: lat, lon, distance along profile (km), depth, stacked amplitude, number of seismic rays intersected

The .dat file containing projected station coordinate, like 'sta_project_newsta.dat'

【3】Finally run 'python wfstretch.py', you need to specify whether you're using RFs with water(crustal multiple) removed and whether you're using new station list or old station list in Line15-16

Output: CCP figure like 'CCP_rwater_newsta.png'



## VII. 'Average' RF (RF_P/radiationstack/)

stack polarity-corrected traces and then calculate 'average' RF



firstly go to ***../***, run 'python rtzcal_event.py -Ssta'

parallel mode: use line 23 instead of line 22 in ***RFcal_multi.py***

output: sac files of rotated components divided by event in [out_path]/[rotate_dir]_event/



then get back to ***RF_P/radiationstack/***

### 1. Direct stack from all traces

run 'python directalign.py', please specify the frequency band in Line 31

output: directstack[freqmin]_[freqmax]hz\_[align(none|maxamp)].png



### 2. Indirect stack

Firstly stack traces of each event respectively, then stack the stacked traces based on different alignment method

run 'python stackRF.py', please specify the frequency band in Line 31

output: test[freqmin]_[freqmax]hz\_[align(none|maxamp)].png



## VIII. SRF (RF_S/)

basically same as the operation in processing PRF

you may adapt the distance range in ***paraSRF.cfg*** for including different phases

while culling the SRFs, we're interested in T channel, so you should use

python PCApickrf.py -Sstation -CT  paraSRF.cfg

