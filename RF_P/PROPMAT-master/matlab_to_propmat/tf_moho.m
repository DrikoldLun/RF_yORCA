function tf = tf_moho(Mohodep,Vs_crust,rayp,varargin)
% TFtype: 1-Moho 2-crustal multiples
TFtype = 1;
samprate = 10;
synthperiod = 0.5;

Pwin = [10 30]; %beforeP, afterP
freq = [0.1 2]; %freqmax, freqmin
RFwin = [10 30]; %beforeP, afterP

narginchk(3,inf);
iv = 1;
while iv+3 < nargin
    switch varargin{iv}
        case 'TFtype'
            TFtype = varargin{iv+1};
        case 'samprate'
            samprate = varargin{iv+1};
        case 'synthperiod'
            synthperiod = varargin{iv+1};
        case 'Pwin'
            Pwin = varargin{iv+1};
        case 'RFwin'
            RFwin = varargin{iv+1};
        case 'freq'
            freq = varargin{iv+1};
    end
    iv = iv+2;
end

model_noMoho = labmodel(0,0,'Vs',[Vs_crust 4.4 4.1]);
model_Moho = labmodel(Mohodep,0,'Vs',[Vs_crust 4.4 4.1]);

[~,~,RF_noMoho] = propmat_syn('rayp',rayp,'zlayt',model_noMoho.zlayt,'zlayb', ...
    model_noMoho.zlayb,'Vp',model_noMoho.Vp,'Vs',model_noMoho.Vs,'rho',model_noMoho.rho, ...
    'samprate',samprate,'synthperiod',synthperiod,'Pwin',Pwin,'RFwin',RFwin,'freq',freq);
[~,~,RF_Moho] = propmat_syn('rayp',rayp,'zlayt',model_Moho.zlayt,'zlayb', ...
    model_Moho.zlayb,'Vp',model_Moho.Vp,'Vs',model_Moho.Vs,'rho',model_Moho.rho, ...
    'samprate',samprate,'synthperiod',synthperiod,'Pwin',Pwin,'RFwin',RFwin,'freq',freq);

RF_time = -RFwin(1):1/samprate:RFwin(2);

if TFtype == 1
    tf = tf_cal(RF_noMoho,RF_Moho,0);
elseif TFtype == 2
    RF_rmul = RF_Moho;
    RF_rmul(find(RF_time>1.5)) = 0;
    tf = tf_cal(RF_rmul,RF_Moho,0);
else
    error('TFtype could only be 1-Moho or 2-crustal multiples!');
end

end