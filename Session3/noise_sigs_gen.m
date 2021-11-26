%% Create speech L
load('../sim_environment/Computed_RIRs.mat')
load('../HRTF.mat') 
source_filename{1} = 'speech1.wav';
speechfilename = {'Speech_Signals/speech1.wav','Speech_Signals/speech2.wav','Speech_Signals/speech3.wav','Speech_Signals/speech4.wav','Speech_Signals/speech5.wav','Speech_Signals/speech6.wav','Speech_Signals/speech7.wav','Speech_Signals/speech8.wav'};
% Number of loudspeakers
J = size(RIR_sources,3);
siglength = 10;
Lh = 400; % Length of impulse responses of listening room
Lg = 2*(Lh-1)/(J-2);      % Length of filters g_j

RIR = RIR_sources(1:Lh,:,:);
% Define the lengths of RIRs and g filters
% Calculate delay for SOE
Delta=ceil(sqrt(room_dim(1)^2+room_dim(2)^2)*fs_RIR/340);

% Define the Toeplitz matrices for left and right ear (left side of SOE)
[~,mic,speaker] = size(RIR_sources);
HL =[];
HR = [];
for j = 1:speaker
    temp_mat = toeplitz(RIR(:,1,j),zeros(Lg,1));
    HL = [HL temp_mat];
end
for j = 1:speaker
    temp_mat = toeplitz(RIR(:,2,j),zeros(Lg,1));
    HR = [HR temp_mat];
end
        

    

% Define the HRTFs for left and right ear (right side of SOE) from the
% loaded HRTF
xL_undelayed = HRTF(:,1); % Left ear
xR_undelayed = HRTF(:,2); % Right ear
xL_delayed = cat(1,zeros(Delta,1),xL_undelayed(Delta+1:Lh,1));
xR_delayed = cat(1,zeros(Delta,1),xR_undelayed(Delta+1:Lh,1));

% Construct H (from HL and HR) and x (from xL and xR) and remove all-zero rows in H, and the corresponding elements in x
H = cat(1,HL,HR);
x = cat(1,xL_delayed,xR_delayed);
[rows,cols]=size(H);
zerows = zeros(1,rows);

idx_nonzeroslines = sum(abs(H),2)> 0;
H_1 = H(idx_nonzeroslines,:);
x_1 = x(idx_nonzeroslines,:);
g = H_1\x_1;

% importing speech signal
speech1 = audioread('../Speech_Signals/speech1.wav');
speech1_res = resample(speech1,8000,44100);
speech1_cut = speech1_res(1:10*8000);

%with standard conv
transf_tot = H_1*g;
transf_L = transf_tot(1:length(transf_tot)/2,:);
transf_R = transf_tot(length(transf_tot)/2+1:end,:);
convL = conv(speech1_cut,transf_L);
convR = conv(speech1_cut,transf_R);
binaural_sig = [convL convR];
%% generate babble noise
noisefilename = {'../Speech_Signals/Babble_noise1.wav','../Speech_Signals/White_noise1.wav'};
noise_amount = size(RIR_noise,3);
[response_length,mic_amount,speaker_amount] = size(RIR_sources);
resample_noise_signals = [];
 %filtered time window in seconds
max_length=siglength*fs_RIR;
for i = 1:noise_amount
    if (i==1) 
        [y_noise,Fs_noise] = audioread(noisefilename{i}); 
        resample_noise = resample(y_noise,fs_RIR,Fs_noise);
        resample_noise_signals = resample_noise(1:max_length,:);
    else
        [y_noise,Fs_noise] = audioread(noisefilename{i}); 
        resample_noise = resample(y_noise,fs_RIR,Fs_noise);
        resample_noise_signals = [resample_noise_signals resample_noise(1:max_length,1)];
    end
end
back_noise = zeros(max_length,mic_amount);
for j = 1:mic_amount
    for i=1:noise_amount
        back_noise(:,j)=...
            back_noise(:,j)+...
            fftfilt(RIR_noise(:,j,i),resample_noise_signals(:,i));
    end
end
%% SNR babble noise
SNR_babbel = 10*log10(var(binaural_sig(:,1))/var(1.15*back_noise(:,1)));
disp(SNR_babbel)
babble_noise_scaled = 1.15*back_noise;
%% generate uncorr noise
uncorr_noise = randn(siglength*Fs_noise,1);
resample_noise = resample(uncorr_noise,fs_RIR,Fs_noise);
resample_noise_signals = resample_noise(1:max_length,:);
back_noise = zeros(max_length,mic_amount);
for j = 1:mic_amount
    for i=1:noise_amount
        back_noise(:,j)=...
            back_noise(:,j)+...
            fftfilt(RIR_noise(:,j,i),resample_noise_signals(:,i));
    end
end
%% SNR uncorr noise
SNR_uncorr = 10*log10(var(binaural_sig(:,1))/var(0.005*back_noise(:,1)));
disp(SNR_uncorr)
uncorr_noise_scaled = 0.005*back_noise;
%% Combine speech sig
combined_sig = binaural_sig(1:80000,:)+uncorr_noise_scaled+babble_noise_scaled;
%% Applying WOLA combined
nfft = 512;
noverlap = 2;
analwin = sqrt(hann(nfft,'periodic'));
synthwin = sqrt(hann(nfft,'periodic')); 
channels = 5;
sizH = size(H_1,2);
L = floor((length(x) - nfft + (nfft / noverlap)) / (nfft / noverlap));
g = -1;
[y_STFT1,f] = WOLA_analysis(combined_sig(:,1),fs_RIR,analwin,nfft,noverlap,g);
[n_STFT1,~] = WOLA_analysis(uncorr_noise_scaled(:,1)+babble_noise_scaled(:,1),fs_RIR,analwin,nfft,noverlap,g);
[x_STFT1,~] = WOLA_analysis(binaural_sig(:,1),fs_RIR,analwin,nfft,noverlap,g);
[y_STFT2,f] = WOLA_analysis(combined_sig(:,2),fs_RIR,analwin,nfft,noverlap,g);
[n_STFT2,~] = WOLA_analysis(uncorr_noise_scaled(:,2)+babble_noise_scaled(:,1),fs_RIR,analwin,nfft,noverlap,g);
[x_STFT2,~] = WOLA_analysis(binaural_sig(:,2),fs_RIR,analwin,nfft,noverlap,g);
y_STFT = cat(3,y_STFT1,y_STFT2);
n_STFT = cat(3,n_STFT1,n_STFT2);
x_STFT = cat(3,x_STFT1,x_STFT2);
%% Plotting spectrogram combined
figure
imagesc(abs(y_STFT1))
ylabel('Freq')
xlabel('Time')
title('Spectrogram combined sig')
colorbar
%% SPP
[noisePowMat, SPP] = spp_calc(combined_sig(:,1),nfft,nfft/noverlap);
[N_freqs, N_frames] = size(y_STFT1);
figure; imagesc(1:N_frames, f,SPP); colorbar; axis xy; set(gcf,'color','w');  
set(gca,'Fontsize',14), xlabel('Time Frames'), ylabel('Frequency (Hz)'), title('Speech Presence Probability for ref mic');
%% MWF
num_mics = 2;

Rnn = cell(N_freqs,1);  Rnn(:) = {1e-6*ones(num_mics,num_mics)};      % Noise Only (NO) corr. matrix. Initialize to small random values
Ryy = cell(N_freqs,1);  Ryy(:) = {1e-6*ones(num_mics,num_mics)};      % Speech + Noise (SPN) corr. matrix. Initialize to small random values
lambda = 0.995;                                                       % Forgetting factors for correlation matrices - can change
SPP_thr = 0.95;                                                       % Threshold for SPP - can change



% For MWF filter for left ear
S_mvdr_mwfL_stft = zeros(N_freqs,N_frames);         
X_mvdr_mwfL_stft = zeros(N_freqs,N_frames);         
N_mvdr_mwfL_stft = zeros(N_freqs,N_frames);  
W_mvdr_mwfL = (1/num_mics)*ones(num_mics,N_freqs); 




% STFT Processing
% Looping through each time frame and each frequency bin
tic
for l=2:N_frames % Time index
    
    for k = 1:N_freqs % Freq index
        
        
        % Create a vector of mic. signals
        Y_kl = squeeze(y_STFT(k,l,1:num_mics));  % M x 1 noisy mic sigs for this freq and time frame
        X_kl = squeeze(x_STFT(k,l,1:num_mics));
        N_kl = squeeze(n_STFT(k,l,1:num_mics));

        % ## Update the correlation matrices using the forgetting factor.
        % Threshold the SPP in order to distinguish between periods of speech and non-speech activity
        speech_tresh = 0.9;
        if SPP(k,l) > speech_tresh
            Ryy{k} = (lambda^2)*Ryy{k}+(1-lambda^2)*(Y_kl*Y_kl');
        else
            Rnn{k} = (lambda^2)*Rnn{k}+(1-lambda^2)*(N_kl*N_kl');
        end
        
        %Only check left mic for treshold, process both mics together
                
                

        
        % Computing the MWF filter using a generalised eigenvalue
        % decomposition of the correlation matrices.
        
        [Q,D] = eig(Ryy{k},Rnn{k});
        s_diff =  diag(D);
        s_diff(1) = 1-(1/s_diff(1));
        s_diff(2:end) =0;
        F = inv(Q')*diag(s_diff)*Q';
          
        
        

       
        W_mvdr_mwfL(:,k) = F(1,:); %Taking the left microphone as a reference
        
        % Filtering the noisy speech, the speech-only, and the noise-only.
        S_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* Y_kl(1:num_mics);
        X_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* X_kl(1:num_mics);
        N_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* N_kl(1:num_mics);
         

        
    end % end freqs
end % end time frames
toc

%% Observe processed STFTs
fs = fs_RIR;
time = 0:1/fs:((length(x)-1)/fs);
clow = -60; chigh = 10;
figure; imagesc(time,f/1000,mag2db(abs(y_STFT(:,:,1))), [clow, chigh]); colorbar; axis xy; set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('microphne signal, 1st mic');
figure; imagesc(time,f/1000,mag2db(abs(S_mvdr_mwfL_stft(:,:,1))), [clow, chigh]); colorbar; axis xy; set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'),title('Enhanced Signal L - MWF');


