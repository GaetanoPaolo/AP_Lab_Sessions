[speech, Fs] =audioread('../Speech_Signals/speech1.wav');
speech = speech(1:5*Fs);
nfft = 2048;
noverlap = 2;
analwin = sqrt(hann(nfft,'periodic'));
synthwin = sqrt(hann(nfft,'periodic')); 
channels = 5;

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

figure(2)
transf = stft(speech, 'Window',analwin,'Overlaplength',D,'FFTLength',nfft);
image(20*log10(abs(transf)))
colormap jet

g = load('g.mat');
g = g.g;
speeches = repmat(speech,1,channels);

% WOLA analysis
[X,f] = WOLA_analysis(speeches,Fs,analwin,nfft,noverlap,g);

%plotting spectrogram
%magn_sq = sum(X,3).*conj(sum(X,3));
magn_sq = sum(X.*conj(X),3);
figure(3)

colormap(hot)
imagesc(magn_sq)
xlabel('Time windows')
ylabel('Freq Bins')
title('WOLA analysis spectrogram')

load('../sim_environment/Computed_RIRs.mat')

%WOLA synthesis
x = WOLA_synthesis(X,synthwin,nfft,noverlap);
%x_res = reshape(x',[1,size(x,1)*size(x,2)]); Is this necessary?
x = repmat(x,1,2);
x_left = x(:,1:channels);
x_right = x(:,channels+1:channels*2);
synth_yeet = [x_left(:,1) x_right(:,1)];
for j = 1:channels
    x_left(:,j) = OLA(x_left(:,j),RIR_sources(1:400,1,j),nfft);
end
for j = 1:channels
    x_right(:,j) = OLA(x_right(:,j),RIR_sources(1:400,2,j),nfft);
end
synth_x = [sum(x_left,2) sum(x_right,2)];


% hold on
% plot(1:length(x_res),speech(1:length(x_res)))
% plot(1:length(x_res),real(x_res))
synth_error = norm(speech(1:length(x_res))-x_res')

