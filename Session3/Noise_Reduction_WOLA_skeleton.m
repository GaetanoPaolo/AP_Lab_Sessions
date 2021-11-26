% Lab 3 for Digital Audio Signal Processing Lab Sessions
% Session 3: Noise reduction in the STFT domain
% R.Ali, G. Bernardi, J.Schott, A. Bertrand
% 2020
%
% The following is the skeleton code for doing the noise reduction in the
% STFT domain. Complete the missing sections accordingly


%% Exercise 3.1: ## Obtain the noisy microphone signals as outlined in the session document.
%
% SOME CONVENTIONS:
%
% Your clean speech signal should be called "speech" and is a binaural signal such that:
% speech = [speech_L speech_R]  % speech_L and speech_R are the clean binaurally synthesised speech signals from session 1
% 
% Your noises should be stored in one binaural variable called "noise"


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Section of code to complete %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

% Create noisy mic signals in the time domain:
y_TD = speech + noise;  % stacked microphone signals


%
%% Apply WOLA analysis to observe signals in the STFT domain, Apply the SPP.


fs = fs_RIR;    % sampling freq
nfft = 512;    % number of DFT points
window = sqrt(hann(nfft,'periodic')); % analysis window
noverlap = 2;   % factor for overlap. noverlap = 2 corresponds to 50% overlap
time = 0:1/fs:((length(x)-1)/fs);


% ## Apply the WOLA analysis to the noisy mic. signals, the speech, and the noise.

[y_STFT,f] = % To complete
[n_STFT,~] = % To complete
[x_STFT,~] = % To complete


% Observe the STFT
clow = -60; chigh = 10; % lower and upper limits for signal power in spectrogram (can change for different resolutions)
[N_freqs, N_frames] = size(y_STFT(:,:,1));

figure; 
imagesc(time, f/1000, mag2db(abs(x_STFT(:,:,1))), [clow, chigh]); colorbar; 
axis xy; set(gca,'fontsize', 14);
set(gcf,'color','w'); xlabel('Time Frame'); ylabel('Frequency (kHz)')




% ## Compute the Speech Presence Probability on the reference microphone
% (you can try the speech-only signal or the noisy-speech in one of the microphones)
% Use the attached spp_calc.m function

[noisePowMat, SPP] =  % To complete


% Observe the SPP
figure; imagesc(1:N_frames, f,SPP); colorbar; axis xy; set(gcf,'color','w');  
set(gca,'Fontsize',14), xlabel('Time Frames'), ylabel('Frequency (Hz)'), title('Speech Presence Probability for ref mic');



%%  Exercise 3.2: ## Implementation of the MWF
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Section of code to complete (3 lines) %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        % Computing the MWF filter using a generalised eigenvalue
        % decomposition of the correlation matrices.
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Section of code to complete ~ 10 lines %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

       
        W_mvdr_mwfL(:,k) =  % Final expression for filter
        
        % Filtering the noisy speech, the speech-only, and the noise-only.
        S_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* Y_kl(1:num_mics);
        X_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* X_kl(1:num_mics);
        N_mvdr_mwfL_stft(k,l) = W_mvdr_mwfL(:,k)'* N_kl(1:num_mics);
        

        
    end % end freqs
end % end time frames
toc





%% Observe processed STFTs

figure; imagesc(time,f/1000,mag2db(abs(y_STFT(:,:,1))), [clow, chigh]); colorbar; axis xy; set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'), title('microphne signal, 1st mic');
figure; imagesc(time,f/1000,mag2db(abs(S_mvdr_mwfL_stft(:,:,1))), [clow, chigh]); colorbar; axis xy; set(gcf,'color','w');set(gca,'Fontsize',14); xlabel('Time (s)'), ylabel('Frequency (Hz)'),title('Enhanced Signal L - MWF');


%% Apply the synthesis stage of the WOLA framework to obtain the time domain equivalents:

s_mwfL = % To complete (time-domain version of S_mvdr_mwfL_stft)
x_mwfL = % To complete (time-domain version of X_mvdr_mwfL_stft) 
n_mwfL = % To complete (time-domain version of N_mvdr_mwfL_stft)


% PLOT SIGNALS
% LISTEN TO SIGNALS!



%% EVALUATION

SNR_in = % Compute input SNR
SNR_out = % Compute output SNR


