function [RF_Time, RF] = IDRF(P,D,dt,t_bounds,alpha,accept_mis,itmax,isplot)
% [RF_Time, RF] = IDRF(P,D,dt,t_bounds,gauss_t,accept_mis,itmax)
%
% Iterative Deconvolution and Receiver-Function Estimation in time domain
%
% P for parent phase (S in Sp case), D for daughter phase.
% dt = time delta (1/samprate) (s)
% t_bounds = time bounds to preserve in eventual RF (really just cut out past t_max)
% gauss_t = 1 std for the gaussian convolved with the delta spikes (recommend 1 s)
% accept_mis = acceptable misfit for (D - conv(P,rf))/length(D)
% itmax = max number of iterations
%
% Sent by Karen Fischer, Jan 2019 
% I think developed by Nick Mancinelli in 2017
% adapted by Hannah Krueger and/or Junlin Hua, early 2019
% lightly edited by Z. Eilon, Jan 2019
%close all;
fs = round(1/dt);
%disp(fs);

%{
%windowing
[~,ci] = max(abs(P));
t_last = 80;
if ci-fs*t_last < 1
    P = cat(1,zeros(fs*t_last-ci+1,1),P(1:floor(ci+fs*t_last)));
    D = cat(1,zeros(fs*t_last-ci+1,1),D(1:floor(ci+fs*t_last)));
else
    P = P(floor(ci-fs*t_last):floor(ci+fs*t_last));
    D = D(floor(ci-fs*t_last):floor(ci+fs*t_last));
end

%detrend
P = detrend(P);
D = detrend(D);
%}
%filtering
%{
if(filt)
    P = filtfilt(b,a,P);
    D = filtfilt(b,a,D);
end
%}
t_max = max(t_bounds); %max preserved point in RF\
t_min = min(t_bounds); %min preserved point in RF\

misfit_old = 999999999999999999;
misfit = sqrt(sum((P-D).^2));

RF = zeros(length(P)*2-1,1);

D_cur = D;

itnum = 0;

%corr with P to find spikes
while (itnum <= itmax) && (misfit_old > accept_mis) 
    [amp_corr,t_corr] = xcorr(D_cur,P); %find highest spike
    auto_corr = xcorr(P);
    [~,ind] = max(abs(amp_corr));
    amp_rf = amp_corr(ind)/auto_corr((length(t_corr)+1)/2);
    RF(ind) = RF(ind)+amp_rf;
    D_sub = conv(P,RF,'same');
    D_cur = D - D_sub;
    %plot(D_cur)
    misfit_old = misfit;
    misfit = sqrt(sum(D_cur.^2))/length(D_cur);
    itnum = itnum+1;
end

RF_Time = t_corr*dt;

%RF(RF_Time>t_max)=0;
%RF = RF(RF_Time<=t_max);
%RF_Time = RF_Time(RF_Time<=t_max);

RF(RF_Time>t_max & RF_Time<t_min)=0;
RF = RF(RF_Time<=t_max & RF_Time>=t_min);
RF_Time = RF_Time(RF_Time<=t_max & RF_Time>=t_min);

gauss = gauss_pulse(1/fs,length(RF),alpha);
RF = real(ifft(fft(RF).*gauss));
%{
if gauss_t~=0
    % gauss_len = length(RF);
    % Gauss_win = gausswin(gauss_len,gauss_sig*4);
    gauss_sig = gauss_t/dt;
    x = linspace(-gauss_sig*4,gauss_sig*4,gauss_sig*8);
    Gauss_win = exp(-x.^2/(2*gauss_sig^2));
    RF = conv(RF,Gauss_win,'same');
end
%}
%RF = flipud(RF);
%RF_Time = fliplr(RF_Time);

%plot(dt*t_corr(t_corr<0 & dt*t_corr>-40),RF(t_corr<0 & dt*t_corr>-40))
%plot(RF_Time(RF_Time>-10 & RF_Time<30),RF(RF_Time>-10 & RF_Time<30))
%plot(RF_Time,RF)

if isplot
y1=RF;y2=RF;
thre = 0*max(abs(RF));
y1(RF<=thre) = thre; y2(RF>=-thre) = -thre;
y1b(1:length(y1)) = thre; y2b(1:length(y2)) = -thre;
fig=fill([RF_Time flip(RF_Time)],[y1' flip(y1b)],'r',[RF_Time flip(RF_Time)],[y2' flip(y2b)],'b'); 
set(fig,{'LineStyle'},{'none'}); hold on;
plot(RF_Time,RF,'k','linewidth',1);
grid;
xlabel('t[s]'); %title(['RF',num2str(index)]);
set(gcf,'unit','normalized','position',[0.2,0.2,0.64,0.32]);
%saveas(gcf,['RFfig/',num2str(index),'.png']);
end
end