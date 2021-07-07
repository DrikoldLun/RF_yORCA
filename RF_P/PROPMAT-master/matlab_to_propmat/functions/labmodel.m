function model = labmodel(Mohodep,LABdep,varargin)
% Mohodep,LABdep 0 denotes no that interface
% model must include 3 layers

Vp = [6.5, 8, 7.2];
Vs = [3.7, 4.4, 4.1];
rho = [2.85, 3.4, 3.2];

narginchk(2,inf);
iv = 1;
while iv+2 < nargin
    switch varargin{iv}
        case 'Vp'
            Vp = varargin{iv+1};
        case 'Vs'
            Vs = varargin{iv+1};
        case 'rho'
            rho = varargin{iv+1};
    end
    iv = iv+2;
end

if ~Mohodep
    model.zlayt = [0, LABdep];
    model.zlayb = [LABdep, 100];
    model.Vp = Vp(2:end);
    model.Vs = Vs(2:end);
    model.rho = rho(2:end);
elseif ~LABdep
    model.zlayt = [0, Mohodep];
    model.zlayb = [Mohodep, 100];
    model.Vp = Vp(1:2);
    model.Vs = Vs(1:2);
    model.rho = rho(1:2);
else
    model.zlayt = [0, Mohodep, LABdep];
    model.zlayb = [Mohodep, LABdep, 100];
    model.Vp = Vp;
    model.Vs = Vs;
    model.rho = rho;
end
end