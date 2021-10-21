%% Speech1
[y_speech1, f_speech1] = audioread('../Speech_Signals/speech1.wav');
f_res = 8000;
y_res = resample(y_speech1,f_res,f_speech1);
x = y_res(1:f_res*10);
binaurals_sig1_1 = [x x];
binaurals_sig2_1 = [x 0.5*x];
binaurals_sig3_1 = [x cat(1,zeros(3,1),x(1:end-3,1))];
load('../HRTF.mat');
x_left = fftfilt(HRTF(:,1),x);
x_right = fftfilt(HRTF(:,2),x);
binaurals_sig4_1 = [x_left x_right];
%% Speech2
[y_speech2, f_speech2] = audioread('../Speech_Signals/speech2.wav');
f_res = 8000;
y_res = resample(y_speech2,f_res,f_speech2);
x = y_res(1:f_res*10);
binaurals_sig1_2 = [x x];
binaurals_sig2_2 = [0.5*x x];
binaurals_sig3_2 = [cat(1,zeros(3,1),x(1:end-3,1)) x];
load('../HRTF.mat');
x_left = fftfilt(HRTF(:,1),x);
x_right = fftfilt(HRTF(:,2),x);
binaurals_sig4_2 = [x_right x_left];
%% Copypasta
binaurals_sig1=binaurals_sig1_1+binaurals_sig1_2;
binaurals_sig2=binaurals_sig2_1+binaurals_sig2_2;
binaurals_sig3=binaurals_sig3_1+binaurals_sig3_2;
binaurals_sig4=binaurals_sig4_1+binaurals_sig4_2;
soundsc(binaurals_sig3, f_res)