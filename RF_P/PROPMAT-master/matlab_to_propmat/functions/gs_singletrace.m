function s = gs_singletrace(data,rayp,varargin)
Vp = 6.5;
samprate = 10;
synthperiod = 0.5;
Hlim = [4 8];
Hsearch = [4.5 7.5];
Vslim = [3. 4.2];
Vssearch = [3. 4.];
freq = [0.1 2];
TFtype = 2;
RFwin = [10 30];
Pwin = [10 30];

narginchk(2,inf);
iv = 1;
while iv+2 < nargin
    switch varargin{iv}
        case 'Vp'
            Vp = varargin{iv+1};
        case 'samprate'
            samprate = varargin{iv+1};
        case 'synthperiod'
            synthperiod = varargin{iv+1};
        case 'RFwin'
            Pwin = varargin{iv+1};
        case 'RFwin'
            Pwin = varargin{iv+1};
        case 'freq'
            freq = varargin{iv+1};
        case 'Hlim'
            Hlim = varargin{iv+1};
        case 'Hsearch'
            Hsearch = varargin{iv+1};
        case 'Vslim'
            Vslim = varargin{iv+1};
        case 'Vssearch'
            Vssearch = varargin{iv+1};
        case 'TFtype'
            TFtype = varargin{iv+1};
    end
    iv = iv+2;
end

dH = 0.1; dVs = 0.1;
Hbg = Hlim(1):dH:Hlim(2);
Vsbg = Vslim(1):dVs:Vslim(2);
searchHb = floor((Hsearch(1)-Hbg(1))/dH)+1;
searchHe = floor((Hsearch(2)-Hbg(1))/dH)+1;
searchVsb = floor((Vssearch(1)-Vsbg(1))/dVs)+1;
searchVse = floor((Vssearch(2)-Vsbg(1))/dVs)+1;
s = ones(length(Vsbg),length(Hbg))*10;
RF_time = (-RFwin(1):1/samprate:RFwin(2));

[m, n] = meshgrid(1:length(Vsbg),1:length(Hbg));
[comb(:,1),comb(:,2)] = deal(reshape(m,[],1), reshape(n,[],1));

RF_nomoho = rf_moho(0,0,rayp);

for i = 1:size(comb,1)
    Vsid = comb(i,1); Vs = Vsbg(Vsid);
    Hid = comb(i,2); h = Hbg(Hid);
    disp(sprintf('%f,%f',Vs,h))
    try
        %tf = tf_moho(h,Vs,rayp,'Vp',Vp,'freq',freq,'TFtype',TFtype,'RFwin',RFwin, ...
        %    'samprate',samprate,'synthperiod',synthperiod,'Pwin',Pwin);
        RF_moho = rf_moho(h,Vs,rayp);
        tf = tf_cal(RF_nomoho,RF_moho,0);
        rf_clean = real(ifft(fft(data).*tf));
        %res_rf = data-rf_clean;
        %disp(Vs,h,sqrt(norm(res_rf(find(RFtime>1.5)))/length(find(RFtime>1.5))));
        %disp(length(rf_clean))
        %disp(length(RF_time))
        %rf_clean = rf_clean(find(RF_time>1.5&RF_time<3));
        %disp(rf_clean)
        s(Vsid,Hid) = norm(rf_clean(find(RF_time>1.5&RF_time<3)));
        %s(Vsid,Hid) = sqrt(norm(rf_clean(find(6>RF_time>1.5)))/length(find(6>RF_time>1.5)));
        disp(s(Vsid,Hid))
    catch err
        disp(sprintf('%f,%f err',Vs,h))
    end
end

s_search = s(searchVsb:searchVse,searchHb:searchHe);
min_s = min(s_search(:));
[i_search,j_search] = find(s_search == min_s);
if length(i_search) > 1
    i = i_search(1)+searchVsb-1; j = j_search(1)+searchHb-1;
else
    i = i_search+searchVsb-1; j = j_search+searchHb-1;
end

Vsfinal = Vsbg(i); Hfinal = Hbg(j);
disp(['H=',num2str(Hfinal),', Vs=',num2str(Vsfinal)]);

figure(1)
gca = pcolor(Hbg,Vsbg,s);
set(gca,'LineStyle','none');
load('cyan.mat');
colorbar; colormap(cyan);
rectangle('Position',[Hsearch(1),Vssearch(1),diff(Hsearch),diff(Vssearch)], ...
    'LineWidth',1,'LineStyle','--');

figure(2)
tf = tf_moho(h,Vs,rayp,'Vp',Vp,'freq',freq,'TFtype',TFtype,'RFwin',RFwin, ...
        'samprate',samprate,'synthperiod',synthperiod,'Pwin',Pwin);
rf_clean = real(ifft(fft(data).*tf));
plot(RF_time,data)
hold on
plot(RF_time,rf_clean)
xlim([-1 10])
end