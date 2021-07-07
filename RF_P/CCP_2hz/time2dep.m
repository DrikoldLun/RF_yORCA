function [depth, x_s, x_p] = time2dep(rayp, time)
addpath('Tinycodes');
%% Model information
%Velocity1D='../IASP91.vel';
%VelocityModel = load(Velocity1D,'-ascii');
velmod = load('IASP91.vel');
VelocityModel = velmod;
rayp = skm2srad(rayp); %ray parameter in s/rad

% Depths
%--------------------------------------------------------------------------
Depths = VelocityModel(:,1);
YAxisRange = [0:0.1:100];
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

%% 1D ray tracing:
%--------------------------------------------------------------------------
    
x_s_lst = cumsum((dz./R) ./ sqrt((1./(rayp^2.* (R./Vs).^-2)) - 1));%Pds piercing distance from station in rad
raylength_s = (dz.*R) ./  (sqrt(((R./Vs).^2) - (rayp^2)).* Vs);
x_p_lst = cumsum((dz./R) ./ sqrt((1./(rayp^2.* (R./Vp).^-2)) - 1));%P piercing distance from station in rad
raylength_p = (dz.*R) ./  (sqrt(((R./Vp).^2) - (rayp^2)).* Vp);

% Calculate Pds travel time
%----------------------------------------------------------------------
Tpds = cumsum((sqrt((R./Vs).^2 - rayp^2) - sqrt((R./Vp).^2 - rayp^2)) .* (dz./R));%the travel time versus depth in columns for all RFs
%Tpds = cumsum(dz./R.*(1./((Vs./R).*sqrt(1-(rayp*Vs./R).^2))-(1./((Vp./R).*sqrt(1-(rayp*Vp./R).^2)))));
%disp(Tpds)

depth = interp1(Tpds,YAxisRange,time);
x_s = interp1(Tpds,x_s_lst,time);
x_p = interp1(Tpds,x_p_lst,time);

return