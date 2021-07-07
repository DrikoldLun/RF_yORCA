clc;clear;fclose all;close all;
%% Datapath and information
addpath('Tinycodes');
Stalist='latlon_newsta.dat';
%Datapath='RFdat_0.1-1hz'; 
Datapath='../../RF_data/CCP/RF50_0.10to2.00/RFdat_culled'
velmod = load('IASP91.vel');
Lineloca = [-133.55 -8;
            -132.45 -4];
Network = 'XX';
sampling = 0.1; shift = 10;
dep_step = 0.5;
YAxisRange=(0:dep_step:150);%the objective depth range ????????????
[dis azi] = distance(Lineloca(1,2),Lineloca(1,1),Lineloca(2,2),Lineloca(2,1));%calculate the azi from the start point to the end point

%order = 3;%bandpass filter parameters ??????????????????
%ny = 1/(2*sampling);
%f1 = 0.03;%low frequency limit??????????
%f2 = 0.2;%high frequency limit
%[B,A] = butter(order, [f1 f2]/ny);

[Stations Stalon Stalat]=textread(Stalist,'%s %f %f %*s %*s %*s',-1);
stanumber=length(Stations);
m=1;
while 1
    if m==stanumber+1
    break,end;
disp([num2str(m) '----' 'station: ' char(Stations(m,:)) ' latitude: ' num2str(Stalat(m)) ' longitude: ' num2str(Stalon(m))]);
m=m+1;
end


%% calculate the Ps piercing points for the events at each station
filename=['pierce0.1-2hz_culled.dat'];xst=exist(filename,'file');fid=fopen(filename,'w+');
RFdepth = struct([]);

rf_time = -shift:sampling:30;
for i =1:stanumber
Station=char(Stations(i,:));stalat=Stalat(i);stalon=Stalon(i);
disp(['the ' num2str(i) 'th' ' station' Station '----------------']);
[event phase evlat evlon evdep dis bazi rayp magnitude f0]=textread(fullfile(Datapath,Station,[Network '.' Station 'finallist.dat']),'%s %s %f %f %f %f %f %f %f %f',-1);
EV_num=length(evdep);
rayp = skm2srad(rayp); %ray parameter in s/rad
Piercelat = zeros(length(YAxisRange),EV_num);%????????
Piercelon = zeros(length(YAxisRange),EV_num);%????????
%read RF data
for m = 1:EV_num
filename=event{m};
%disp(char(filename));
datafile=fullfile(Datapath,Station,[filename,'_',Station,'_PRF.R']);
%datar(:,m)=load(datafile);
datar(:,m)= rmparent(load(datafile),rf_time);
end

% bandpass filter ????????:
% for j = 1:EV_num
% datar(:,j) = filtfilt(B,A,datar(:,j));
% end

% 1D ray tracing:
%--------------------------------------------------------------------------
[rfdepth, EndIndex, x_s, x_p] = PSRF2depth(datar, rayp, YAxisRange, sampling, shift, velmod);

% Ps piercing point calculation
%calculate piercing point lat/lon at 70 km, 410 km, 660 km????????????
for l  = 1:EV_num
    %[Pierce_loc(l,1) Pierce_loc(l,2)] =latlon_from(stalat,stalon,bazi(l),deg2km(rad2deg(x_s(821,l))));%70
    %[Pierce_loc(l,3), Pierce_loc(l,4)] =latlon_from(stalat,stalon,bazi(l),deg2km(rad2deg(x_s(1041,l))));%410
    %[Pierce_loc(l,5) Pierce_loc(l,6)] =latlon_from(stalat,stalon,bazi(l),deg2km(rad2deg(x_s(1321,l))));%660
    [Pierce_loc(l,7) Pierce_loc(l,8)] =latlon_from(stalat,stalon,bazi(l),deg2km(rad2deg(x_s(91,l))));
    %[Pierce_loc(l,1) Pierce_loc(l,2)] =latlon_from(stalat,stalon,bazi(l),deg2km(rad2deg(x_s(floor(70/dep_step)+1,l))));%70
    %[Pierce_loc(l,3), Pierce_loc(l,4)] =latlon_from(stalat,stalon,bazi(l),deg2km(rad2deg(x_s(floor(410/dep_step)+1,l))));%410
    %[Pierce_loc(l,5) Pierce_loc(l,6)] =latlon_from(stalat,stalon,bazi(l),deg2km(rad2deg(x_s(floor(660/dep_step)+1,l))));%660
    [Pierce_loc(l,7) Pierce_loc(l,8)] =latlon_from(stalat,stalon,bazi(l),deg2km(rad2deg(x_s(floor(45/dep_step)+1,l))));
end

for l = 1:EV_num    
        [Piercelat(:,l) Piercelon(:,l)] = latlon_from(stalat,stalon,bazi(l),deg2km(rad2deg(x_s(:,l))));
end
%project point to profile for depth
         newlat = [];
         newlon = [];
        
      [dis_center, azi_center] = distance(Lineloca(1,2),Lineloca(1,1),Piercelat,Piercelon);
        newlat=[];newlon=[];
        dis_along = atand(tand(dis_center)).*cosd(azi - azi_center); 
        for m=1:size(dis_along,2)
        [newlat(:,m),newlon(:,m)] =latlon_from(Lineloca(1,2),Lineloca(1,1),azi,deg2km(dis_along(:,m)));
        end


%save to the structure
RFdepth(i).Station = Station;
RFdepth(i).stalat = stalat;
RFdepth(i).stalon = stalon;
RFdepth(i).Depthrange = YAxisRange;
RFdepth(i).events = event;
RFdepth(i).phases = phase;
RFdepth(i).moveout_correct = rfdepth;
RFdepth(i).Piercelat = Piercelat;
RFdepth(i).Piercelon = Piercelon;
RFdepth(i).projectlat = newlat; %Piercelat;
RFdepth(i).projectlon = newlon; %Piercelon;
RFdepth(i).StopIndex = EndIndex;
RFdepth(i).Piercepoint = [Pierce_loc(:,1),Pierce_loc(:,2),Pierce_loc(:,3),Pierce_loc(:,4),Pierce_loc(:,5), Pierce_loc(:,6)];
%output piercing locations to a file
fprintf(fid,'%s\n',Station);

for m =1:EV_num
fprintf(fid,'%s %f %f %f %f %f %f %f %f\n',char(event(m,:)),Pierce_loc(m,1),Pierce_loc(m,2),Pierce_loc(m,3),Pierce_loc(m,4),Pierce_loc(m,5), Pierce_loc(m,6),Pierce_loc(m,7), Pierce_loc(m,8));
end

clear datar Pierce_loc Piercelat Piercelon;
end
save RFdepth_culled_newsta.mat RFdepth
fclose(fid);

function rfdata_new = rmparent(rfdata,rftime)
rfdata_afterP = rfdata(find(rftime>=0));
rftime_afterP = rftime(find(rftime>=0));
parentend = 0;
for i = 2:length(rfdata_afterP)-1
    if (rfdata_afterP(i-1) >= rfdata_afterP(i)) && (rfdata_afterP(i+1) >= rfdata_afterP(i))
        parentend = rftime_afterP(i);
        break
    end
end
rfdata_new = rfdata(1:end);
rfdata_new(find(rftime<parentend)) = 0;
end