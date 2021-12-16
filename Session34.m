%% Session 3 code
%% Single channel left ear
addpath("../Session2")
load('../sim_environment/Computed_RIRs_session3.mat')
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
% SNR uncorr noise
SNR_uncorr = 10*log10(var(binaural_sig(:,1))/var(0.005*back_noise(:,1)));
uncorr_noise_scaled = 0.005*back_noise;
% Combine speech sig
combined_sig = binaural_sig(1:80000,:)+uncorr_noise_scaled+babble_noise_scaled;
% Applying WOLA combined
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
% Single channel left ear
% SPP
[noisePowMat, SPP] = spp_calc(combined_sig(:,1),nfft,nfft/noverlap);%binaural_sig(1:80000,1)
[N_freqs, N_frames] = size(y_STFT1);
% MWF
num_mics = 1; %When using single channel, change this to 1, result is very poor in quality (-16dB SNR)

Rnn = cell(N_freqs,1);  Rnn(:) = {(10e-3)*ones(num_mics,num_mics)};      % Noise Only (NO) corr. matrix. Initialize to small random values
Ryy = cell(N_freqs,1);  Ryy(:) = {(10e-3)*ones(num_mics,num_mics)};      % Speech + Noise (SPN) corr. matrix. Initialize to small random values

lambda_n = 0.6;
lambda_y = 0.7;                                                       % Forgetting factors for correlation matrices - can change
SPP_thr = 0.99;                                                       % Threshold for SPP - can change




% For MWF filter for left ear
S_mvdr_mwfL_stft = zeros(N_freqs,N_frames);         
X_mvdr_mwfL_stft = zeros(N_freqs,N_frames);         
N_mvdr_mwfL_stft = zeros(N_freqs,N_frames);  
W_mvdr_mwfL = (1/num_mics)*ones(num_mics,N_freqs); 

count = sum(SPP > SPP_thr);
noisy_frames = sum(count==0);
SN_mvdr_mwfL_stft = zeros(N_freqs,noisy_frames);         
n_ind = 2;
SS_mvdr_mwfL_stft = zeros(N_freqs,N_frames-noisy_frames);
s_ind = 2;

% STFT Processing
for l=2:N_frames % Time index
    
    for k = 1:N_freqs % Freq index
        
        
        % Create a vector of mic. signals
        Y_kl = squeeze(y_STFT(k,l,1:num_mics));  % M x 1 noisy mic sigs for this freq and time frame
        X_kl = squeeze(x_STFT(k,l,1:num_mics));
        N_kl = squeeze(n_STFT(k,l,1:num_mics));

        % ## Update the correlation matrices using the forgetting factor.
        % Threshold the SPP in order to distinguish between periods of speech and non-speech activity
        
        if SPP(k,l) > SPP_thr
            Ryy{k} = (lambda_y^2)*Ryy{k}+(1-lambda_y^2)*(Y_kl*Y_kl');
        else
            Rnn{k} = (lambda_n^2)*Rnn{k}+(1-lambda_n^2)*(Y_kl*Y_kl');
        end
        
        %Only check left mic for treshold, process both mics together
                

        
        % Computing the MWF filter using a generalised eigenvalue
        % decomposition of the correlation matrices.
        
        [V,D] = eig(Ryy{k},Rnn{k});
        [s_diff,I] =  sort(diag(D),'descend');
        Q = inv(V)';
        Q = Q(:,I);
        s_diff(1) = 1-(1/s_diff(1));
        s_diff(2:end) =0;
        F = (Q')\(diag(s_diff)*Q');
          
        W_mvdr_mwfL(:,k) = F(:,1); %Taking the left microphone as a reference

        % Filtering the noisy speech, the speech-only, and the noise-only.
        S_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* Y_kl(1:num_mics);
        X_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* X_kl(1:num_mics);
        N_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* N_kl(1:num_mics);
        if SPP(k,l) < SPP_thr
            SN_mvdr_mwfL_stft(k,n_ind) = W_mvdr_mwfL(:,k)'*Y_kl(1:num_mics);
        else
            SS_mvdr_mwfL_stft(k,s_ind) = W_mvdr_mwfL(:,k)'*Y_kl(1:num_mics);
        end 
        
    end % end freqs
    if(count(l)==0)
        n_ind = n_ind + 1;
    else
        s_ind = s_ind + 1;
    end
end % end time frames

% Observe processed STFTs
fs = fs_RIR;
time = 0:1/fs:((length(x)-1)/fs);
clow = -60; chigh = 10;
figure; 
subplot(2,1,1)
imagesc(time,f/1000,mag2db(abs(y_STFT(:,:,1))), [clow, chigh]); colorbar; axis xy; set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('Microphone signal, Left ear Single channel');
subplot(2,1,2)
imagesc(time,f/1000,mag2db(abs(S_mvdr_mwfL_stft(:,:,1))), [clow, chigh]); colorbar; axis xy; set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'),title('Enhanced Signal L - MWF Single channel');

% WOLA synthesis
y_filt_left = WOLA_synthesis(S_mvdr_mwfL_stft,analwin,nfft,noverlap);
yn_filt_left = WOLA_synthesis(SN_mvdr_mwfL_stft,analwin,nfft,noverlap);
yy_filt_left = WOLA_synthesis(SS_mvdr_mwfL_stft,analwin,nfft,noverlap);
n_filt_left = WOLA_synthesis(N_mvdr_mwfL_stft,analwin,nfft,noverlap);
x_filt_left = WOLA_synthesis(X_mvdr_mwfL_stft,analwin,nfft,noverlap);
% SNR computation
noise_sum = uncorr_noise_scaled(:,1)+babble_noise_scaled(:,1);
SNR_in = mag2db(rssq(binaural_sig(:,1))/rssq(noise_sum(:,1)));
disp('SNR input left ear single channel MWF')
disp(SNR_in)
SNR_out = mag2db(rssq(yy_filt_left(:))/rssq(yn_filt_left(:)));
disp('SNR output left ear single channel MWF')
disp(SNR_out)

%% Single channel right ear
% SPP
[noisePowMat, SPP] = spp_calc(combined_sig(:,2),nfft,nfft/noverlap);%binaural_sig(1:80000,2)
[N_freqs, N_frames] = size(y_STFT2);
% MWF
num_mics = 1; %When using single channel, change this to 1, result is very poor in quality (-16dB SNR)

Rnn = cell(N_freqs,1);  Rnn(:) = {(10e-3)*ones(num_mics,num_mics)};      % Noise Only (NO) corr. matrix. Initialize to small random values
Ryy = cell(N_freqs,1);  Ryy(:) = {(10e-3)*ones(num_mics,num_mics)};      % Speech + Noise (SPN) corr. matrix. Initialize to small random values

lambda_n = 0.6;
lambda_y = 0.7;                                                       % Forgetting factors for correlation matrices - can change
SPP_thr = 0.99;                                                       % Threshold for SPP - can change




% For MWF filter for left ear
S_mvdr_mwfR_stft = zeros(N_freqs,N_frames);         
X_mvdr_mwfR_stft = zeros(N_freqs,N_frames);         
N_mvdr_mwfR_stft = zeros(N_freqs,N_frames);  
W_mvdr_mwfR = (1/num_mics)*ones(num_mics,N_freqs); 

count = sum(SPP > SPP_thr);
noisy_frames = sum(count==0);
SN_mvdr_mwfL_stft = zeros(N_freqs,noisy_frames);         
n_ind = 2;
SS_mvdr_mwfL_stft = zeros(N_freqs,N_frames-noisy_frames);
s_ind = 2;

% STFT Processing
% Looping through each time frame and each frequency bin
for l=2:N_frames % Time index
    
    for k = 1:N_freqs % Freq index
        
        
        % Create a vector of mic. signals
        Y_kl = squeeze(y_STFT(k,l,2));  % M x 1 noisy mic sigs for this freq and time frame
        X_kl = squeeze(x_STFT(k,l,2));
        N_kl = squeeze(n_STFT(k,l,2));

        % ## Update the correlation matrices using the forgetting factor.
        % Threshold the SPP in order to distinguish between periods of speech and non-speech activity
        
        if SPP(k,l) > SPP_thr
            Ryy{k} = (lambda_y^2)*Ryy{k}+(1-lambda_y^2)*(Y_kl*Y_kl');
        else
            Rnn{k} = (lambda_n^2)*Rnn{k}+(1-lambda_n^2)*(Y_kl*Y_kl');
        end
        
        %Only check left mic for treshold, process both mics together
                

        
        % Computing the MWF filter using a generalised eigenvalue
        % decomposition of the correlation matrices.
        
        [V,D] = eig(Ryy{k},Rnn{k});
        [s_diff,I] =  sort(diag(D),'descend');
        Q = inv(V)';
        Q = Q(:,I);
        s_diff(1) = 1-(1/s_diff(1));
        s_diff(2:end) =0;
        F = (Q')\(diag(s_diff)*Q');
          
        W_mvdr_mwfR(:,k) = F(:,1); %Taking the left microphone as a reference

        % Filtering the noisy speech, the speech-only, and the noise-only.
        S_mvdr_mwfR_stft(k,l) = W_mvdr_mwfR(:,k)'* Y_kl(1:num_mics);
        X_mvdr_mwfR_stft(k,l) = W_mvdr_mwfR(:,k)'* X_kl(1:num_mics);
        N_mvdr_mwfR_stft(k,l) = W_mvdr_mwfR(:,k)'* N_kl(1:num_mics);
        if SPP(k,l) < SPP_thr
            SN_mvdr_mwfR_stft(k,n_ind) = W_mvdr_mwfR(:,k)'*Y_kl(1:num_mics);
        else
            SS_mvdr_mwfR_stft(k,s_ind) = W_mvdr_mwfR(:,k)'*Y_kl(1:num_mics);
        end 
        
    end % end freqs
    if(count(l)==0)
        n_ind = n_ind + 1;
    else
        s_ind = s_ind + 1;
    end
end % end time frames

% Observe processed STFTs
fs = fs_RIR;
time = 0:1/fs:((length(x)-1)/fs);
clow = -60; chigh = 10;
figure; 
subplot(2,1,1)
imagesc(time,f/1000,mag2db(abs(y_STFT(:,:,2))), [clow, chigh]); colorbar; axis xy; set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('Microphone signal, Right ear Single channel');
subplot(2,1,2)
imagesc(time,f/1000,mag2db(abs(S_mvdr_mwfR_stft(:,:,1))), [clow, chigh]); colorbar; axis xy; set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'),title('Enhanced Signal R - MWF Single channel');

% WOLA synthesis
y_filt_right = WOLA_synthesis(S_mvdr_mwfR_stft,analwin,nfft,noverlap);
yn_filt_right = WOLA_synthesis(SN_mvdr_mwfR_stft,analwin,nfft,noverlap);
yy_filt_right = WOLA_synthesis(SS_mvdr_mwfR_stft,analwin,nfft,noverlap);
n_filt_right = WOLA_synthesis(N_mvdr_mwfR_stft,analwin,nfft,noverlap);
x_filt_right = WOLA_synthesis(X_mvdr_mwfR_stft,analwin,nfft,noverlap);
% SNR computation
noise_sum = uncorr_noise_scaled(:,2)+babble_noise_scaled(:,2);
SNR_in = mag2db(rssq(binaural_sig(:,2))/rssq(noise_sum(:,1)));
disp('SNR input right ear single channel MWF')
disp(SNR_in)
SNR_out = mag2db(rssq(yy_filt_right(:))/rssq(yn_filt_right(:)));
disp('SNR output right ear single channel MWF')
disp(SNR_out)

%% Multi channel left & right ears
% SPP
[noisePowMat, SPP] = spp_calc(combined_sig(:,1),nfft,nfft/noverlap);%binaural_sig(1:80000,1)
[N_freqs, N_frames] = size(y_STFT1);

% MWF
num_mics = 2; %When using single channel, change this to 1, result is very poor in quality (-16dB SNR)

Rnn = cell(N_freqs,1);  Rnn(:) = {(10e-3)*ones(num_mics,num_mics)};      % Noise Only (NO) corr. matrix. Initialize to small random values
Ryy = cell(N_freqs,1);  Ryy(:) = {(10e-3)*ones(num_mics,num_mics)};      % Speech + Noise (SPN) corr. matrix. Initialize to small random values

lambda_n = 0.6;
lambda_y = 0.7;                                                       % Forgetting factors for correlation matrices - can change
SPP_thr = 0.9;                                                       % Threshold for SPP - can change




% For MWF filter for left ear
S_mvdr_mwfL_stft = zeros(N_freqs,N_frames);         
X_mvdr_mwfL_stft = zeros(N_freqs,N_frames);         
N_mvdr_mwfL_stft = zeros(N_freqs,N_frames);  
W_mvdr_mwfL = (1/num_mics)*ones(num_mics,N_freqs); 
W_mvdr_mwfR = (1/num_mics)*ones(num_mics,N_freqs); 
S_mvdr_mwfR_stft = zeros(N_freqs,N_frames); 

count = sum(SPP > SPP_thr);
noisy_frames = sum(count==0);
SN_mvdr_mwfL_stft = zeros(N_freqs,noisy_frames);         
n_ind = 2;
SS_mvdr_mwfL_stft = zeros(N_freqs,N_frames-noisy_frames);
s_ind = 2;

% STFT Processing
% Looping through each time frame and each frequency bin
for l=2:N_frames % Time index
    
    for k = 1:N_freqs % Freq index
        
        
        % Create a vector of mic. signals
        Y_kl = squeeze(y_STFT(k,l,1:num_mics));  % M x 1 noisy mic sigs for this freq and time frame
        X_kl = squeeze(x_STFT(k,l,1:num_mics));
        N_kl = squeeze(n_STFT(k,l,1:num_mics));

        % ## Update the correlation matrices using the forgetting factor.
        % Threshold the SPP in order to distinguish between periods of speech and non-speech activity
        
        if SPP(k,l) > SPP_thr
            Ryy{k} = (lambda_y^2)*Ryy{k}+(1-lambda_y^2)*(Y_kl*Y_kl');
        else
            Rnn{k} = (lambda_n^2)*Rnn{k}+(1-lambda_n^2)*(Y_kl*Y_kl');
        end
        
        %Only check left mic for treshold, process both mics together
                

        
        % Computing the MWF filter using a generalised eigenvalue
        % decomposition of the correlation matrices.
        
        [V,D] = eig(Ryy{k},Rnn{k});
        [s_diff,I] =  sort(diag(D),'descend');
        Q = inv(V)';
        Q = Q(:,I);
        s_diff(1) = 1-(1/s_diff(1));
        s_diff(2:end) =0;
        F = (Q')\(diag(s_diff)*Q');

       
        W_mvdr_mwfL(:,k) = F(:,1); %Taking the left microphone as a reference
        W_mvdr_mwfR(:,k) = F(:,2);

        % Filtering the noisy speech, the speech-only, and the noise-only.
        S_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* Y_kl(1:num_mics);
        X_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* X_kl(1:num_mics);
        N_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* N_kl(1:num_mics);
        if SPP(k,l) < SPP_thr
            SN_mvdr_mwfL_stft(k,n_ind) = W_mvdr_mwfL(:,k)'*Y_kl(1:num_mics);
        else
            SS_mvdr_mwfL_stft(k,s_ind) = W_mvdr_mwfL(:,k)'*Y_kl(1:num_mics);
        end

        S_mvdr_mwfR_stft(k,l) = W_mvdr_mwfR(:,k)'*Y_kl(1:num_mics); 

    end % end freqs
    if(count(l)==0)
        n_ind = n_ind + 1;
    else
        s_ind = s_ind + 1;
    end
end % end time frames

% Observe processed STFTs
fs = fs_RIR;
time = 0:1/fs:((length(x)-1)/fs);
clow = -60; chigh = 10;
figure; 
subplot(2,1,1)
imagesc(time,f/1000,mag2db(abs(y_STFT(:,:,1))), [clow, chigh]); colorbar; axis xy; set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('Microphone signal, Left ear Multi-channel');
subplot(2,1,2)
imagesc(time,f/1000,mag2db(abs(S_mvdr_mwfL_stft(:,:,1))), [clow, chigh]); colorbar; axis xy; set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'),title('Enhanced Signal L - MWF');

figure; 
subplot(2,1,1)
imagesc(time,f/1000,mag2db(abs(y_STFT(:,:,2))), [clow, chigh]); colorbar; axis xy; set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('Microphone signal, Right ear Multi-channel');
subplot(2,1,2)
imagesc(time,f/1000,mag2db(abs(S_mvdr_mwfR_stft(:,:,1))), [clow, chigh]); colorbar; axis xy; set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'),title('Enhanced Signal R - MWF');

% WOLA synthesis
y_filt_left = WOLA_synthesis(S_mvdr_mwfL_stft,analwin,nfft,noverlap);
y_filt_right = WOLA_synthesis(S_mvdr_mwfL_stft,analwin,nfft,noverlap);
yn_filt_left = WOLA_synthesis(SN_mvdr_mwfL_stft,analwin,nfft,noverlap);
yy_filt_left = WOLA_synthesis(SS_mvdr_mwfL_stft,analwin,nfft,noverlap);
n_filt_left = WOLA_synthesis(N_mvdr_mwfL_stft,analwin,nfft,noverlap);
x_filt_left = WOLA_synthesis(X_mvdr_mwfL_stft,analwin,nfft,noverlap);
% SNR computation
noise_sum = uncorr_noise_scaled(:,1)+babble_noise_scaled(:,1);
SNR_in = mag2db(rssq(binaural_sig(:,1))/rssq(noise_sum(:,1)));
disp('SNR input left ear multi channel MWF')
disp(SNR_in)
SNR_out = mag2db(rssq(yy_filt_left(:))/rssq(yn_filt_left(:)));
disp('SNR output left ear multi channel MWF')
disp(SNR_out)


%% Session 4 ANC: FxLMS Noise only
% Fx-NLMS 
clear all;
% close all

addpath("../Session2")

% Load RIRs
load('../sim_environment/Computed_RIRs_session4.mat');
% Set length
sigLenSec = 10;

speakersel = 3;


%
% Plot the RIRs of the noise
[len,~] = size(RIR_noise);
[len2,~,~] = size(RIR_sources);
% Generate noisy mic signal

% Read in the noise source (resample if necessary)
[y_noise,Fs_noise] = audioread('../Speech_Signals/White_noise1.wav'); 
resample_noise = resample(y_noise,fs_RIR,Fs_noise);
filt_noise = fftfilt(RIR_noise(:,1),resample_noise);
filt_noise = filt_noise(1:sigLenSec*fs_RIR);

% Generate desired speech 
% importing speech signal
speech1 = audioread('../Speech_Signals/speech1.wav');
speech1_res = resample(speech1,8000,44100);
filtspeech_res = fftfilt(RIR_sources(:,1,speakersel),speech1_res);
filtspeech_res = filtspeech_res(1:sigLenSec*fs_RIR);
speech_cut = speech1_res(1:sigLenSec*fs_RIR);


% FxLMS Noise only


M = 80;  % Length of secondary path (RIR)
L = 600;  % Adaptive filter length

mu = 0.14;   % Step size
delta = 5*10^(-5);

W = zeros(L,1); % Initialize adaptive filter

sigLenSample = sigLenSec*fs_RIR;
speech_cut = cat(1,zeros(L+M-1,1),speech_cut);
filtspeech_res = cat(1,zeros(L+M-1,1),filtspeech_res);
x = cat(1,zeros(L+M-1,1),resample_noise(1:sigLenSample));
e = zeros(1,sigLenSample);
h = RIR_sources(1:M,1,speakersel);
d = filt_noise;

SNR_input = 10*log10(var(filtspeech_res)/var(filt_noise));
SNR_inputbis = 10*log10(var(speech_cut)/var(x));

for n = 1:sigLenSample

    % STEP 1 : Arrange the previous L + M − 1 samples of x(n) up until the
    % current sample, n, in an [M x L] Hankel matrix X_Hmat (this will be
    % use for the filtering operation)
    xseg = x(n:n+L+M-1); % current processed filtered noise segment
    temp = flip(xseg);
    half1 = temp(1:L)';
    half2 = temp(L:end-1)';
    X_Hmat = hankel(half1,half2);
    % STEP 2: Filter the signal x(n) with the adaptive filter w(n) to
    % generate loudspeaker signals y(n). Store the samples of y(n) into the
    % vector y


    y =  W'*X_Hmat; 
    % STEP 3: Compute the error signal e(n)
    e(n) = d(n) + h'*y'; 

    % STEP 4: Pre-filter x(n) with the (estimated) RIRs to obtain the
    % filtered x(n), \bar{x}(n). Store the samples of \bar{x}(n) into
    % the vector xf

    xf = X_Hmat*h;
    %xf = fftfilt(h,xseg);
    % STEP 5: Update the filter w(n). Store the samples of w(n) into
    % the vector w

    W = W - (mu/(norm(xf)^2+delta))*xf*e(n);
end

%% Plotting Fx-NLMS outputs
figure
title("Left ear noise only")
hold on 
plot(1:sigLenSample,filt_noise)
plot(1:sigLenSample,e)
legend('d(n)','e(n)')
hold off


%% FxLMS Noise + Audio


M = 80;  % Length of secondary path (RIR)
L = 600;  % Adaptive filter length

mu = 0.14;   % Step size
delta = 5*10^(-5);

W = zeros(L,1); % Initialize adaptive filter

sigLenSample = sigLenSec*fs_RIR;
speech_cut = cat(1,zeros(L+M-1,1),speech_cut);
filtspeech_res = cat(1,zeros(L+M-1,1),filtspeech_res);
x = cat(1,zeros(L+M-1,1),resample_noise(1:sigLenSample));
e = zeros(1,sigLenSample);
h = RIR_sources(1:M,1,speakersel);
d = filt_noise;


for n = 1:sigLenSample

    % STEP 1 : Arrange the previous L + M − 1 samples of x(n) up until the
    % current sample, n, in an [M x L] Hankel matrix X_Hmat (this will be
    % use for the filtering operation)
    xseg = x(n:n+L+M-1); % current processed filtered noise segment
    temp = flip(xseg);
    half1 = temp(1:L)';
    half2 = temp(L:end-1)';
    X_Hmat = hankel(half1,half2);
    % STEP 2: Filter the signal x(n) with the adaptive filter w(n) to
    % generate loudspeaker signals y(n). Store the samples of y(n) into the
    % vector y


    y =  W'*X_Hmat; 
    % STEP 3: Compute the error signal e(n)
    e(n) = d(n) + h'*y'+ filtspeech_res(n); 
    
    % STEP 4: Pre-filter x(n) with the (estimated) RIRs to obtain the
    % filtered x(n), \bar{x}(n). Store the samples of \bar{x}(n) into
    % the vector xf

    xf = X_Hmat*h;
    %xf = fftfilt(h,xseg);
    % STEP 5: Update the filter w(n). Store the samples of w(n) into
    % the vector w

    W = W - (mu/(norm(xf)^2+delta))*xf*e(n);
end

% Plotting Fx-NLMS outputs
figure
title("Left ear noise + audio")
hold on 
plot(1:sigLenSample,filt_noise)
plot(1:sigLenSample,e)
hold off
legend('d(n)','e(n)')


%% MFx-NLMS Noise Only

% Load RIRs
load('../sim_environment/Computed_RIRs_session4circ.mat');
% Set length
sigLenSec =10;
% Number of speakers
speakers = size(RIR_sources,3);
AUDIO_SPEECH_3D = 0;
%
[len,J] = size(RIR_noise);

% Generate noisy mic signal

% Read in the noise source (resample if necessary)
[y_noise,Fs_noise] = audioread('../Speech_Signals/White_noise1.wav'); 
resample_noise = resample(y_noise,fs_RIR,Fs_noise);
filt_noiseL = fftfilt(RIR_noise(:,1),resample_noise);
filt_noiseL = filt_noiseL(1:sigLenSec*fs_RIR);
filt_noiseR = fftfilt(RIR_noise(:,2),resample_noise);
filt_noiseR = filt_noiseR(1:sigLenSec*fs_RIR);
filt_noise = [filt_noiseL filt_noiseR];

% Generate binaural sig
if AUDIO_SPEECH_3D
    % Load measured HRTFs
    load('../HRTF.mat') 

    % Define the signal length
    siglength = 10;


    % Define the lengths of RIRs and g filters
    Lh = 901; % Length of impulse responses of listening room
    Lg = 2*(Lh-1)/(speakers-2);      % Length of filters g_j

    % Truncate impulse response to reduce computational complexity
    RIR = RIR_sources(1:Lh,:,:);

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

    % Solve the SOE
    g = H_1\x_1;
    % Generate desired speech 
    % importing speech signal
    speech1 = audioread('../Speech_Signals/speech1.wav');
    speech1_res = resample(speech1,8000,44100);
    %filtspeech_res = fftfilt(RIR_sources(:,1,speakersel),speech1_res);
    %filtspeech_res = filtspeech_res(1:sigLenSec*fs_RIR);
    speech_cut = speech1_res(1:sigLenSec*fs_RIR);

    transf_tot = H_1*g;
    transf_L = transf_tot(1:length(transf_tot)/2,:);
    transf_R = transf_tot(length(transf_tot)/2+1:end,:);
    convL = fftfilt(transf_L,speech_cut);
    convR = fftfilt(transf_R,speech_cut);
    binaural_sig = [convL convR];
    % Achieve 0 dB between left binaural sig and  left noise
    SNR_left = 1;
    count = 1;
    while SNR_left > 0.001
        count  = count +1;
        SNR_left = snr((0.9999^(count))*convL, filt_noiseL);
    end
    sig_scale = 0.9999^count;
    binaural_sig_scaled = [convL*sig_scale, convR*sig_scale];
end
% MFxLMS

M = 80;  % Length of secondary path (RIR)
L = 600;  % Adaptive filter length

%mu = 0.5;   % Step size
mu = 0.1;
delta = 5*10^(-5);

sigLenSample = sigLenSec*fs_RIR;
x = cat(1,zeros(L+M-1,1),resample_noise(1:sigLenSample));
e_L = zeros(1,sigLenSample);
e_R = zeros(1,sigLenSample);
hL = zeros(M,speakers);
hR = zeros(M,speakers);
for i = 1:speakers
    hL(:,i) = RIR_sources(1:M,1,i);
    hR(:,i) = RIR_sources(1:M,2,i);
end

W = zeros(L,speakers);
d = filt_noise;

tic
for n=1:sigLenSample
    
    % STEP 1 : Arrange the previous L + M − 1 samples of x(n) up until the
    % current sample, n, in an [M x L] Hankel matrix X_Hmat (this will be
    % use for the filtering operation). Store the samples in the
    % appropriate form.
    xseg = x(n:n+L+M-1); % current processed filtered noise segment
    temp = flip(xseg);
    half1 = temp(1:L)';
    half2 = temp(L:end-1)';
    X_Hmat = hankel(half1,half2);
    

    % STEP 2: Filter the signal x(n) with the adaptive filter w(n) to
    % generate loudspeaker signals y(n). Store the samples in
    % the appropriate form.
    
    y = zeros(speakers,M);
    for k = 1:speakers
        y(k,:) = W(:,k)'*X_Hmat;
    end
    
    % STEP 3: Compute the error signals e_L(n) and e_R(n). Store them in
    % the appropriate form.

    e_L(n) = d(n,1);
    e_R(n) = d(n,2);
    for k = 1:speakers
        e_L(n) = e_L(n)+hL(:,k)'*y(k,:)';
        e_R(n) = e_R(n)+hR(:,k)'*y(k,:)';
    end
    if AUDIO_SPEECH_3D
        e_L(n) = e_L(n)+ binaural_sig_scaled(n,1);
        e_R(n) = e_R(n)+ binaural_sig_scaled(n,2);
    end

    % STEP 4: Pre-filter x(n) with the (estimated) RIRs to obtain the
    % filtered x(n), \bar{x}(n). Store the samples in the appropriate form.
    
    xf = zeros(L,speakers,J);
    for s = 1:speakers
        xf(:,s,1) = X_Hmat*hL(:,s);
        xf(:,s,2) = X_Hmat*hR(:,s);
    end

    % STEP 5: Update the filter w(n). Store the samples in the appropriate
    % form.
    for s = 1:speakers
        E = [e_L(n) e_R(n)];
        W(:,s) = W(:,s) - (mu/(norm(squeeze(xf(:,s,:)),"fro")^2+delta))*squeeze(xf(:,s,:))*E';
    end
end
toc

% Plotting Fx-NLMS outputs
figure
subplot(2,1,1)
hold on
plot(1:sigLenSample,filt_noise(:,1))
plot(1:sigLenSample,e_L)
if AUDIO_SPEECH_3D
    plot(1:sigLenSample, binaural_sig_scaled(:,1))
end
hold off
if AUDIO_SPEECH_3D
    legend('d(n)','e(n)','binaural left')
else
    legend('d(n)','e(n)')
end
title('Left filtered noise and error')
xlabel('Amplitude')
ylabel('Time')
subplot(2,1,2)
hold on 
plot(1:sigLenSample,filt_noise(:,2))
plot(1:sigLenSample,e_R)
if AUDIO_SPEECH_3D
    plot(1:sigLenSample, binaural_sig_scaled(:,2))
end
hold off
if AUDIO_SPEECH_3D
    legend('d(n)','e(n)','binaural right')
else
    legend('d(n)','e(n)')
end
title('Right filtered noise and error')
xlabel('Amplitude')
ylabel('Time')


%% MFx-NLMS Noise + speech
% Load RIRs
load('../sim_environment/Computed_RIRs_session4circ.mat');
% Set length
sigLenSec =10;
% Number of speakers
speakers = size(RIR_sources,3);
AUDIO_SPEECH_3D = 1;
%%
[len,J] = size(RIR_noise);

%% Generate noisy mic signal

% Read in the noise source (resample if necessary)
[y_noise,Fs_noise] = audioread('../Speech_Signals/White_noise1.wav'); 
resample_noise = resample(y_noise,fs_RIR,Fs_noise);
filt_noiseL = fftfilt(RIR_noise(:,1),resample_noise);
filt_noiseL = filt_noiseL(1:sigLenSec*fs_RIR);
filt_noiseR = fftfilt(RIR_noise(:,2),resample_noise);
filt_noiseR = filt_noiseR(1:sigLenSec*fs_RIR);
filt_noise = [filt_noiseL filt_noiseR];

%% Generate binaural sig
if AUDIO_SPEECH_3D
    % Load measured HRTFs
    load('../HRTF.mat') 

    % Define the signal length
    siglength = 10;


    % Define the lengths of RIRs and g filters
    Lh = 901; % Length of impulse responses of listening room
    Lg = 2*(Lh-1)/(speakers-2);      % Length of filters g_j

    % Truncate impulse response to reduce computational complexity
    RIR = RIR_sources(1:Lh,:,:);

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

    % Solve the SOE
    g = H_1\x_1;
    % Generate desired speech 
    % importing speech signal
    speech1 = audioread('../Speech_Signals/speech1.wav');
    speech1_res = resample(speech1,8000,44100);
    %filtspeech_res = fftfilt(RIR_sources(:,1,speakersel),speech1_res);
    %filtspeech_res = filtspeech_res(1:sigLenSec*fs_RIR);
    speech_cut = speech1_res(1:sigLenSec*fs_RIR);

    transf_tot = H_1*g;
    transf_L = transf_tot(1:length(transf_tot)/2,:);
    transf_R = transf_tot(length(transf_tot)/2+1:end,:);
    convL = fftfilt(transf_L,speech_cut);
    convR = fftfilt(transf_R,speech_cut);
    binaural_sig = [convL convR];
    % Achieve 0 dB between left binaural sig and  left noise
    SNR_left = 1;
    count = 1;
    while SNR_left > 0.001
        count  = count +1;
        SNR_left = snr((0.9999^(count))*convL, filt_noiseL);
    end
    sig_scale = 0.9999^count;
    binaural_sig_scaled = [convL*sig_scale, convR*sig_scale];
end
% MFxLMS

M = 80;  % Length of secondary path (RIR)
L = 600;  % Adaptive filter length

%mu = 0.5;   % Step size
mu = 0.1;
delta = 5*10^(-5);

sigLenSample = sigLenSec*fs_RIR;
x = cat(1,zeros(L+M-1,1),resample_noise(1:sigLenSample));
e_L = zeros(1,sigLenSample);
e_R = zeros(1,sigLenSample);
hL = zeros(M,speakers);
hR = zeros(M,speakers);
for i = 1:speakers
    hL(:,i) = RIR_sources(1:M,1,i);
    hR(:,i) = RIR_sources(1:M,2,i);
end

W = zeros(L,speakers);
d = filt_noise;

tic
for n=1:sigLenSample
    
    % STEP 1 : Arrange the previous L + M − 1 samples of x(n) up until the
    % current sample, n, in an [M x L] Hankel matrix X_Hmat (this will be
    % use for the filtering operation). Store the samples in the
    % appropriate form.
    xseg = x(n:n+L+M-1); % current processed filtered noise segment
    temp = flip(xseg);
    half1 = temp(1:L)';
    half2 = temp(L:end-1)';
    X_Hmat = hankel(half1,half2);
    

    % STEP 2: Filter the signal x(n) with the adaptive filter w(n) to
    % generate loudspeaker signals y(n). Store the samples in
    % the appropriate form.
    
    y = zeros(speakers,M);
    for k = 1:speakers
        y(k,:) = W(:,k)'*X_Hmat;
    end
    
    % STEP 3: Compute the error signals e_L(n) and e_R(n). Store them in
    % the appropriate form.

    e_L(n) = d(n,1);
    e_R(n) = d(n,2);
    for k = 1:speakers
        e_L(n) = e_L(n)+hL(:,k)'*y(k,:)';
        e_R(n) = e_R(n)+hR(:,k)'*y(k,:)';
    end
    if AUDIO_SPEECH_3D
        e_L(n) = e_L(n)+ binaural_sig_scaled(n,1);
        e_R(n) = e_R(n)+ binaural_sig_scaled(n,2);
    end

    % STEP 4: Pre-filter x(n) with the (estimated) RIRs to obtain the
    % filtered x(n), \bar{x}(n). Store the samples in the appropriate form.
    
    xf = zeros(L,speakers,J);
    for s = 1:speakers
        xf(:,s,1) = X_Hmat*hL(:,s);
        xf(:,s,2) = X_Hmat*hR(:,s);
    end

    % STEP 5: Update the filter w(n). Store the samples in the appropriate
    % form.
    for s = 1:speakers
        E = [e_L(n) e_R(n)];
        W(:,s) = W(:,s) - (mu/(norm(squeeze(xf(:,s,:)),"fro")^2+delta))*squeeze(xf(:,s,:))*E';
    end
end
toc

% Plotting Fx-NLMS outputs
figure
subplot(2,1,1)
hold on
plot(1:sigLenSample,filt_noise(:,1))
plot(1:sigLenSample,e_L)
if AUDIO_SPEECH_3D
    plot(1:sigLenSample, binaural_sig_scaled(:,1))
end
hold off
if AUDIO_SPEECH_3D
    legend('d(n)','e(n)','binaural left')
else
    legend('d(n)','e(n)')
end
title('Left filtered speech + noise and error')
xlabel('Amplitude')
ylabel('Time')
subplot(2,1,2)
hold on 
plot(1:sigLenSample,filt_noise(:,2))
plot(1:sigLenSample,e_R)
if AUDIO_SPEECH_3D
    plot(1:sigLenSample, binaural_sig_scaled(:,2))
end
hold off
if AUDIO_SPEECH_3D
    legend('d(n)','e(n)','binaural right')
else
    legend('d(n)','e(n)')
end
title('Right filtered speech + noise and error')
xlabel('Amplitude')
ylabel('Time')

% How does the convergence of the algorithm compare to the noise only case? How do you
%expect this would be affected in practice?: In noise only time sequences,
%the convergence is similar to the noise only case. However, it is visible
%that the added 3D speech has an amplifying effect when the speech is
%acitve, especially towards the end of an utterance (see
%where_it_goes_bad). This might result in distortion/cracking noises at the
%end of words in practice. The cause might be the strong errors at the
%beginning of the word which have repercussions on the weights and cause
%regular overshoots of the noise + binaural speech at the end of the
%word.Lowering the learning rate mu from 0.5 to 0.1 strongly reduces this
%effect (as the abrupt changes at start of the speech are smoothed), but
%this results in a slower convergence at the beginnning of the audio.







