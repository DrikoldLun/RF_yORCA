function [R_pm,Z_pm,RF_pm] = propmat_syn(varargin)
%% Make synthetic receiver function, plot
%% init
model = struct('nlay',2,...
                      'zlayt',[0 5]',...
                      'zlayb',[  5 30.]',...
                      'Vs',   [3.5 4.1]');
model.Vp = [6. 7.2]';
model.rho = [2.72 3.03]';

samprate = 10;
rayp = 0.06;
synthperiod = 0.5;

isplot = 0;

Pwin = [10 30]; %beforeP, afterP
freq = [0.1 2]; %freqmax, freqmin
RFwin = [10 30]; %beforeP, afterP

ID = [];

narginchk(0,inf);
iv = 1;
while iv < nargin
    switch varargin{iv}
        case 'zlayt'
            model.zlayt = varargin{iv+1}'; %row vector
        case 'zlayb'
            model.zlayb = varargin{iv+1}';
        case 'Vp'
            model.Vp = varargin{iv+1}';
        case 'Vs'
            model.Vs = varargin{iv+1}';
        case 'rho'
            model.rho = varargin{iv+1}';
        case 'samprate'
            samprate = varargin{iv+1};
        case 'rayp'
            rayp = varargin{iv+1};
        case 'synthperiod'
            synthperiod = varargin{iv+1};
        case 'isplot'
            isplot = varargin{iv+1};
        case 'Pwin'
            Pwin = varargin{iv+1};
        case 'RFwin'
            RFwin = varargin{iv+1};
        case 'freq'
            freq = varargin{iv+1};
        case 'ID'
            ID = [varargin{iv+1},''];
    end
    iv = iv+2;
end

model.nlay = length(model.zlayt);
%% CA structure - match Zhu and Kanamori 2000

%rayp = sind(inc)./model.Vp(end);
inc = asind(rayp*model.Vp(end));
[trP,ttP] = run_propmat(model,ID,'Ps',samprate,inc,synthperiod);
% output is TRZ, switch to ZRT, normalise, flip to positive
dataP = trP(:,[3,1,2])./maxab(maxab(trP));

%% compute arrival of parent, shift traces to this zero
t_arP =  ttP(dataP(:,1)==max(abs(dataP(:,1))));
ttP = ttP - t_arP;

%% check output
chans = {'Z','R','T'};

%result = py.tr.Ps_t_amp(modCA_ZK2000.Vp,modCA_ZK2000.Vs,modCA_ZK2000.zlayb-modCA_ZK2000.zlayt,modCA_ZK2000.rho,deg2rad(inc),0);

T1 = 1/freq(2);
T2 = 1/freq(1);
[b_pm,a_pm]=butter(2,[1/samprate/T2,1/samprate/T1]);

R_pm = dataP(:,2);
Z_pm = dataP(:,1);
R_pm = filtfilt(b_pm,a_pm,R_pm);
Z_pm = filtfilt(b_pm,a_pm,Z_pm);


amp_pm = max(Z_pm);
band_info = sprintf('%.2f to %.2f[Hz]',1/T2,1/T1);

R_pm = cutP(R_pm,amp_pm,samprate,Pwin(1),Pwin(2));
Z_pm = cutP(Z_pm,amp_pm,samprate,Pwin(1),Pwin(2));

alpha = 5;
RF_info = sprintf('alpha=%d',alpha);

[RF_time_pm,RF_pm] = IDRF(Z_pm,R_pm,1/samprate,[-RFwin(1),RFwin(2)],alpha,1e-4,1e3,false);

if isplot
%plotZR(RF_pm,1,[1,1],1,samprate,2,10,'red','');
figure(1)
plot(RF_time_pm,RF_pm,'color','red','linewidth',1.5)
xlim([-2 20])
hold on
end
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
end
