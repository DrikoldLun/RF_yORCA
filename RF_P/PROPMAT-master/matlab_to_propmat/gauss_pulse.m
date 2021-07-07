function gauss = gauss_pulse(dt, nt, alpha)
df = 1./(dt*nt);
nf_mid = ceil(0.5*nt);
w = 2*pi*df*[0:nf_mid-1];

gauss = zeros(nt,1);
gauss(1:nf_mid) = exp(-(0.5*w./alpha).^2)./dt;
gauss(end-nf_mid+1:end) = flipud(gauss(1:nf_mid));
%gauss(1:end) = gauss(1:end)/max(ifft(gauss(1:end)));