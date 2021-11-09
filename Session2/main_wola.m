[speech, Fs] =audioread('../Speech_Signals/speech1.wav');
speech = speech(1:5*Fs);
nfft = 512;
noverlap = 2;
analwin = sqrt(hann(nfft,'periodic'));
synthwin = sqrt(hann(nfft,'periodic')); 

%checking for perfect reconstruction
N = nfft;
d = 2;
D = N/d;
sum_array = zeros(1,D);
for i = 1:D
    sum_win = 0;
    for k = 0:d-1
        sum_win = sum_win + analwin(i+k*D)*synthwin(i+k*D);
    end
    sum_array(i) = sum_win;
end

if sum(sum_array) == D
    disp('Perfect reconstruction')
else
    disp('No perfect reconstruction')
end
figure(1)
spectrogram(speech,analwin,noverlap,nfft)
title('Spectrogram function')
g = load('g.mat');
g = g.g;
speeches = repmat(speech,5);
% WOLA analysis
[X,f] = WOLA_analysis(speeches,Fs,analwin,nfft,noverlap,g);

%plotting spectrogram
% magn_sq = abs(X).^2;
% figure(2)
% colormap(hot)
% imagesc(magn_sq)
% xlabel('Time windows')
% ylabel('Freq Bins')
% title('WOLA analysis spectrogram')

%WOLA synthesis
x = WOLA_synthesis(X,synthwin,nfft,noverlap);
%x_res = reshape(x',[1,size(x,1)*size(x,2)]); Is this necessary?
hold on
plot(1:length(x_res),speech(1:length(x_res)))
plot(1:length(x_res),real(x_res))
synth_error = norm(speech(1:length(x_res))-x_res')

