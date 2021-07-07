addpath('functions')
%{
RF_stack_dir = 'PCAstack/RF50_0.10to2.00';
RF_stack_path = fullfile(RF_stack_dir,['WC03','_stack.dat']);
RF_stack = load(RF_stack_path);
stack_t = RF_stack(:,1);
stack_amp = RF_stack(:,2);
s = gs_singletrace(stack_amp(1:5:end),0.05);
%}
%{
rf = rf_moho(6,3.7,0.06);
rf1 = rf_moho(5.9,3.6,0.06);
rf_time = -10:0.1:30;
tf = tf_moho(6,3.7,0.06);
rf_cleaned = ifft(fft(rf1).*tf);
plot(rf_time,rf1)
hold on
plot(rf_time,rf_cleaned)
xlim([-1 10])
%}
RF_stack_dir = 'PCAstack/RF50_0.10to2.00';
RF_stack_path = fullfile(RF_stack_dir,['WC03','_PCA.dat']);
RF_stack = load(RF_stack_path);
stack_t = RF_stack(:,1);
stack_amp = RF_stack(:,2);
stack_amp = stack_amp(1:5:end);

rf_time = -10:0.1:30;
%tf = tf_moho(5.3,3.5,0.06);
rf_syn = rf_moho(6.3,3.5,0.06);
rf0 = rf_moho(0,0,0.06);
tf = tf_cal(rf0/max(rf0)*max(stack_amp),rf_syn/max(rf_syn)*max(stack_amp),0);
rf_cleaned = real(ifft(fft(stack_amp).*tf));
a1=plot(rf_time,stack_amp); label1 = 'obs';
hold on
a2=plot(rf_time,rf_cleaned); label2 = 'cleaned';
%a3=plot(rf_time,rf_syn); label3 = 'syn';
legend([a1;a2], label1, label2);
xlim([-1 10])