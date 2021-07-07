%% Make synthetic receiver function, plot
clear all

T1 = 0.5;
T2 = 10;
samprate = 25;
[b_pm,a_pm]=butter(2,[1/samprate/T2,1/samprate/T1]);
t_before = 3;
t_after = 15;

R_raw = load('syn_semd/raw.r.semd');
Z_raw = load('syn_semd/raw.z.semd');
R_water = load('syn_semd/water.r.semd');
Z_water = load('syn_semd/water.z.semd');
R_raw = R_raw(:,2);
Z_raw = Z_raw(:,2);
R_water = R_water(:,2);
Z_water = Z_water(:,2);

R_raw = filtfilt(b_pm,a_pm,R_raw);
Z_raw = filtfilt(b_pm,a_pm,Z_raw);
R_water = filtfilt(b_pm,a_pm,R_water);
Z_water = filtfilt(b_pm,a_pm,Z_water);

amp_raw = 1;
amp_water = 1;
R_raw = cutP(R_raw,amp_raw,samprate,5,15);
Z_raw = cutP(Z_raw,amp_raw,samprate,5,15);
R_water = cutP(R_water,amp_water,samprate,5,15);
Z_water = cutP(Z_water,amp_water,samprate,5,15);

band_info = sprintf('%.2f to %.2f[Hz]',1/T2,1/T1);

alpha = 5;
RF_info = sprintf('alpha=%d',alpha);

[RF_time_raw,RF_raw] = IDRF(Z_raw,R_raw,1/samprate,[-t_before,t_after],alpha,1e-4,1e3,false);
[RF_time_water,RF_water] = IDRF(Z_water,R_water,1/samprate,[-t_before,t_after],alpha,1e-4,1e3,false);
%figure(2)
%plot(RF_time_sem,RF_sem,RF_time_pm,RF_pm)
%plotZR(RF_sem,1,'RF',fs_sem,2,10,'red',RF_info);
plotZR(Z_raw,1,[3,1],1,samprate,2,10,'black','Z');
plotZR(Z_water,1,[3,1],1,samprate,2,10,'red','Z');
legend('raw','water')
plotZR(R_raw,1,[3,1],2,samprate,2,10,'black','R');
plotZR(R_water,1,[3,1],2,samprate,2,10,'red','R');
legend('raw','water')
plotZR(RF_raw,1,[3,1],3,samprate,2,10,'black','PRF');
plotZR(RF_water,1,[3,1],3,samprate,2,10,'red','PRF');
legend('raw','water')

%% ---------------------------- SUBFUNCTIONS ----------------------------
%% ---------------------------- SUBFUNCTIONS ----------------------------
%% ---------------------------- SUBFUNCTIONS ----------------------------

function data_tmp = plotZR(data,norm,rowcol,pos,fs,t_before,t_after,color,title_add)
    figure(1)
    set(gca,'xlim',[-t_before t_after])
    set(gcf,'pos',[177 139 800 799]);
    [amp,nc] = max(data);
    n_before = t_before*fs;
    n_after = t_after*fs;
    data_tmp = data(nc-n_before:nc+n_after);
    data_tmp = data_tmp/norm;
    t = [-t_before:1/fs:t_after];
    %{
    if strcmp(channel,'Z')
        subplot(3,1,1); hold on
        plot(t,data_tmp,'color',color)
    end
    if strcmp(channel,'R')
        subplot(3,1,2); hold on
        plot(t,data_tmp,'color',color)
    end
    if strcmp(channel,'RF')
        subplot(3,1,3); hold on
        plot(t,data_tmp,'color',color)
    end
    %}
    %xlim([-t_before t_after])
    subplot(rowcol(1),rowcol(2),pos); hold on
    plot(t,data_tmp,'color',color,'linewidth',1.5)
    title(title_add)
    xlabel('t[s]')
end

function data_tmp = cutP(data,norm,fs,t_before,t_after)
    [amp,nc] = max(data);
    n_before = t_before*fs;
    n_after = t_after*fs;
    data_tmp = data(nc-n_before:nc+n_after);
    data_tmp = data_tmp/norm;
end

function [ B,ind ] = maxab( A )
    % [ B,ind ] = maxab( A )
    %   This function returns the largest value in data series A (or in each
    %   column, if A is a matrix. This will be the largest absolute value,
    %   irrespective if that is negative or positive

    if all(isnan(A))
        B = nan; ind = [];
        return
    end

    if isrow(A), A = A(:); end

    mabA = max(abs(A));
    abA = abs(A);
    B = zeros(1,size(A,2));
    ind = zeros(1,size(A,2));
    for ii = 1:size(A,2)
        B(ii) = unique(A(abA(:,ii)==mabA(ii),ii));
        ind(ii) = find(A(:,ii)==B(ii),1,'first');
    end

end


function rho  = sed_vs2rho( Vs )
% rho  = sed_vs2rho( Vs )
%   Empirical scaling of Vs to rho for sedimentary rocks.
%   Equations from Shen and Ritzwoller (JGR, 2016) equation (2), based
%   on results of Brocher (BSSA 2005).

rho = 1.227 + 1.53*Vs - 0.837*Vs.^2 + 0.207*Vs.^3 - 0.0166*Vs.^4;

end

function rho  = mantle_vs2rho( Vs,Zkm )
% rho  = mantle_vs2rho( Vs,Zkm )
%   Empirical scaling of Vs to rho for mantle rocks. density is scaled from
%   Vs using empirical scaling over a range of upper mantle P,T, where the
%   rho and Vs values are computed for "basal Lherzolite" using the
%   calculator of Hacker and Abers 2016 at conditions betwen 50 and 300 km,
%   with a mantle potential temperature of 1300 C, geotherm of 18C/GPa and
%   P calculated as Z(km)/32. From these values, we compute a best fitting
%   scaling function for Vs/rho as a function of pressure along the
%   geotherm and show that for temperature heterogeneity of �60 C the error
%   in computed rho is less than 0.3%. See empirical_VsRho.m

tf = 1.337 + ((175-Zkm)/125)*0.0141;

rho = Vs./tf;

end

function [ Vsv,Vsh ] = VsvVsh_from_VsXi( Vs,xi )
%[ Vsv,Vsh ] = VsvVsh_from_VsXi( Vs,Xi )
%   Function to calculate Vsv and Vsh from the voigt average velocity (Vs)
%   and the value of xi, which describes the radial anisotropy, where
%
%   Vs^2 = (Vsh^2 + 2Vsv^2)/3
%   xi = Vsh^2/Vsv^2 (=N/L)


Vsv = Vs .* sqrt(3./(xi+2));
Vsh = Vs .* sqrt(3.*xi./(xi + 2));


end
