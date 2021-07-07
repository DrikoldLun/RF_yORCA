clc;clear;fclose all;close all;
addpath('Tinycodes');
load RFdepth_rwater_newsta.mat
%Out_path = 'ccp_0.1-1hz.dat';
Out_path = 'ccp_0.1-0.5hz_rwater_newsta.dat';
Lineloca = [-133.55 -8;
            -132.45 -4];
YAxisRange = RFdepth(1).Depthrange;
Stack_range = (0:0.5:150);%????????
domperiod = 5;%dominant period of Ps waves for PRF
Profile_width = 4;
[dis,azi] = distance(Lineloca(1,2),Lineloca(1,1),Lineloca(2,2),Lineloca(2,1));%calculate the azi from the start point to the end point
Profile_range = (0:Profile_width:deg2km(dis));
[Profile_lat, Profile_lon] = latlon_from(Lineloca(1,2),Lineloca(1,1),azi,Profile_range);
RFlength = length(Stack_range);
Stack_RF = zeros(length(Stack_range),length(Profile_range));
Event_count = zeros(length(Stack_range),length(Profile_range));
%% Model information
Velocity1D='IASP91.vel';
VelocityModel = load(Velocity1D,'-ascii');


% Depths
%--------------------------------------------------------------------------
Depths = VelocityModel(:,1);
% Velocities
%--------------------------------------------------------------------------
Vp = VelocityModel(:,2);
Vs = VelocityModel(:,3);
% Interpolate velocity model to match depth range and increments
%--------------------------------------------------------------------------
Vp = interp1(Depths,Vp,YAxisRange)';
Vs = interp1(Depths,Vs,YAxisRange)';
Depths = YAxisRange';

% Depth intervals
%--------------------------------------------------------------------------
dz = [0; diff(Depths)];
% Radial shells
%--------------------------------------------------------------------------
R = 6371 - Depths;
%%
%find the events located within each circle bin and stack:



% hhh= waitbar(0,'Please wait---depth...');
fid = fopen(Out_path, 'w+');
for i = 1:length(Profile_range)
    disp(['calculate the RF stacks at the distance of ' num2str(Profile_range(i)) ' km along the profile-------'])
    %fprintf(fid, '>\n');
  for j =1:length(Stack_range)
%        disp(['----calculate the RF stacks at the depth of ' num2str(Stack_range(j)) ' km -------'])
    for k =1:length(RFdepth)    
        [row,col] = find(distance(Profile_lat(i),Profile_lon(i),...
            RFdepth(k).projectlat(2*Stack_range(j) + 1,:),RFdepth(k).projectlon(2*Stack_range(j) + 1,:))...
             <= km2deg(sqrt(0.5*domperiod*Vs(2*Stack_range(j) + 1)*Stack_range(j))));
%         <= km2deg(45));%??????????70km?????????????? 
        if ~isempty(col)
            Stack_RF(j,i) = Stack_RF(j,i) + sum(RFdepth(k).moveout_correct(2*Stack_range(j) + 1,col));
            Event_count(j,i) = Event_count(j,i) + length(col);
        end
    end
    if Event_count(j,i) > 0
    Stack_RF(j,i) = Stack_RF(j,i)/Event_count(j,i);    
    end
    fprintf(fid, '%f %f %f %f %f %d\n', Profile_lat(i),Profile_lon(i),Profile_range(i),Stack_range(j),Stack_RF(j,i),Event_count(j,i));
  end
%   Stack_data{i,1} = Profile_lat(i); Stack_data{i,2} = Profile_lon(i); Stack_data{i,3} = Profile_range(i);
%   Stack_data{i,4} = Stack_RF(:,i); Stack_data{i,5} = Event_count(:,i);
end
fclose(fid);
% B=strcat(out_pathname,out_filename);
% save(B,'Stack_data')

fid = fopen('sta_project_newsta.dat', 'w+');
for i = 1:length(RFdepth)
    dis = distance(RFdepth(i).projectlat(1,1),RFdepth(i).projectlon(1,1),Lineloca(1,2),Lineloca(1,1));
    %fprintf(fid, '%s %f %f\n',RFdepth(i).Station,RFdepth(i).projectlon(1,1),RFdepth(i).projectlat(1,1));
    fprintf(fid, '%s %f\n',RFdepth(i).Station,deg2km(dis));
end
fclose(fid);
