%% Make synthetic receiver function, plot
clear all
% close all

%% CA structure - match Zhu and Kanamori 2000
modCA_ZK2000 = struct('nlay',2,...
                      'zlayt',[0 5.2]',...
                      'zlayb',[  5.2 30.]',...
                      'Vs',   [3.5 4.1]');
modCA_ZK2000.Vp = [6. 7.2]';
modCA_ZK2000.rho = [2.72 3.03]';

%% use PropMAT to model Ps arrivals
samprate = 25;
inc = 15;
synthperiod = 1;
rayp = sind(inc)./modCA_ZK2000.Vp(end);
[trP,ttP] = run_propmat(modCA_ZK2000,[],'Ps',samprate,inc,synthperiod);
% output is TRZ, switch to ZRT, normalise, flip to positive
dataP = trP(:,[3,1,2])./maxab(maxab(trP));

%% compute arrival of parent, shift traces to this zero
t_arP =  ttP(dataP(:,1)==max(abs(dataP(:,1))));
ttP = ttP - t_arP;


%% check output
chans = {'Z','R','T'};
figure(44);clf; set(gcf,'pos',[177 139 1478 799]);
for ic = 1:3
    % P
    subplot(3,1,ic)
    plot(ttP,dataP(:,ic),'b','linewidth',2)

    set(gca,'ylim',[-1 1],'xlim',[-30 30],'fontsize',13)
    ylabel(chans{ic},'fontsize',20)
    if ic == 3, xlabel('Time since P arrival (s)','fontsize',20); end
    if ic == 1, title('P-wave (s)','fontsize',25); end

end

result = py.tr.Ps_t_amp(modCA_ZK2000.Vp,modCA_ZK2000.Vs,modCA_ZK2000.zlayb-modCA_ZK2000.zlayt,modCA_ZK2000.rho,deg2rad(inc),0);
T1 = 1;
T2 = 10;
[b_pm,a_pm]=butter(2,[1/samprate/T2,1/samprate/T1]);
t_before = 3;
t_after = 15;

R_sem = load('semd/AA.S0015.BXX.semd');
Z_sem = load('semd/AA.S0015.BXZ.semd');
RF_stack_dir = 'RFstack/RF50_0.10to1.00';
%RF_stack_dir = 'RFstack/RF50_0.10to1.00';
RF_stack_path = fullfile(RF_stack_dir,'CC04_stack.dat');
RF_stack = load(RF_stack_path);
stack_t = RF_stack(:,1);
stack_amp = RF_stack(:,2);
fs_sem = round(1/(R_sem(2,1) - R_sem(1,1)));
[b_sem,a_sem]=butter(2,[1/fs_sem/T2,1/fs_sem/T1]);
R_sem = R_sem(:,2);
Z_sem = Z_sem(:,2);
R_sem = filtfilt(b_sem,a_sem,R_sem);
Z_sem = filtfilt(b_sem,a_sem,Z_sem);
R_pm = dataP(:,2);
Z_pm = dataP(:,1);
R_pm = filtfilt(b_pm,a_pm,R_pm);
Z_pm = filtfilt(b_pm,a_pm,Z_pm);


amp_sem = max(Z_sem);
amp_pm = max(Z_pm);
band_info = sprintf('%.2f to %.2f[Hz]',1/T2,1/T1);
%R_sem = plotZR(R_sem,amp_sem,'R',fs_sem,2,10,'red',band_info);
%Z_sem = plotZR(Z_sem,amp_sem,'Z',fs_sem,2,10,'red',band_info);
R_pm = plotZR(R_pm,amp_pm,'R',samprate,2,10,'blue',band_info);
Z_pm = plotZR(Z_pm,amp_pm,'Z',samprate,2,10,'blue',band_info);
alpha = 5;
RF_info = sprintf('alpha=%d',alpha);
%[RF_time_sem,RF_sem] = IDRF(Z_sem,R_sem,1/fs_sem,[-t_before,t_after],alpha,1e-4,1e3,false);
[RF_time_pm,RF_pm] = IDRF(Z_pm,R_pm,1/samprate,[-t_before,t_after],alpha,1e-4,1e3,false);
%figure(2)
%plot(RF_time_sem,RF_sem,RF_time_pm,RF_pm)
%plotZR(RF_sem,1,'RF',fs_sem,2,10,'red',RF_info);
plotZR(RF_pm,1,'RF',samprate,2,10,'blue',RF_info);
%plotZR(RF_sem,1,'RF',fs_sem,2,10,'red');
%plotZR(RF_pm,1,'RF',samprate,2,10,'blue');
%{
figure(2); hold on
plot(RF_time_sem,RF_sem,'r-',RF_time_pm,RF_pm,'b-')
grid
legend({'SEM' 'Propmat'})
%}

figure(1)

subplot(3,1,3)
hold on
stack_amp = stack_amp(find(-2<=stack_t<=10));
stack_t = stack_t(find(-2<=stack_t<=10));
plot(stack_t,stack_amp,'color','black')
xlim([-2 10])
legend('Syn','Obs')

for i = 1:4
subplot(3,1,1)
hold on
plot(result{'t'}{i}-result{'t'}{1},result{'Zamp'}{i}/result{'Zamp'}{1},'k.')
subplot(3,1,2)
hold on
plot(result{'t'}{i}-result{'t'}{1},result{'Ramp'}{i}/result{'Zamp'}{1},'k.')
subplot(3,1,3)
plot(result{'t'}{i}-result{'t'}{1},0,'ro')
end
return

%% ---------------------------- SUBFUNCTIONS ----------------------------
%% ---------------------------- SUBFUNCTIONS ----------------------------
%% ---------------------------- SUBFUNCTIONS ----------------------------

function data_tmp = plotZR(data,norm,channel,fs,t_before,t_after,color,title_add)
    figure(1)
    set(gca,'xlim',[-t_before t_after],'fontsize',13)
    set(gcf,'pos',[177 139 800 799]);
    [amp,nc] = max(data);
    n_before = t_before*fs;
    n_after = t_after*fs;
    data_tmp = data(nc-n_before:nc+n_after);
    data_tmp = data_tmp/norm;
    t = [-t_before:1/fs:t_after];
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
    %xlim([-t_before t_after])
    title([channel,' (',title_add,')'])
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
