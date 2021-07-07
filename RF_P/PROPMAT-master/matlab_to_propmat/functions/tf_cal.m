function tf = tf_cal(rf_cleaned,rf_original,c)
% freq spectrum of RFs
f_rf_c = fft(rf_cleaned);
f_rf_o = fft(rf_original);
    
%numerator & denominator of transfer function
tf_t = f_rf_c.*conj(f_rf_o);
tf_b = f_rf_o.*conj(f_rf_o);
    
% water level
wl = max(real(tf_b))*c; 
tf_b(find(tf_b<wl)) = wl;
    
% transfer function
tf = tf_t./tf_b;
end