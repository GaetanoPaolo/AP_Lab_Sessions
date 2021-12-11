% Session 4 - Exercise 2 - Multi-channel FxNLMS ANC
%           - Exercise 3 - Multi-channel FxNLMS ANC in 3D audio scenario
% 
%
% Main points:
% (1) Generate the noisy microphone signals
% (2) Implement the ANC.
% (3) Compute input/output SNRs and SNR improvement
% (4) Implement the filter for left and right ears. Hear the synthesised signal


clear all;
% close all

addpath("../Session2")

% Load RIRs
load('../sim_environment/Computed_RIRs_session4.mat');
% Set length
sigLenSec =10;
% Number of speakers
speakers = size(RIR_sources,3);
AUDIO_SPEECH_3D = 1;
%%
% Plot the RIRs of the noise
[len,J] = size(RIR_noise);
figure(1); clf;

subplot(2,1,1)
plot(1:len,RIR_noise(:,1))
subplot(2,1,2)
plot(1:len,RIR_noise(:,2))
%% Generate noisy mic signal

% Read in the noise source (resample if necessary)
[y_noise,Fs_noise] = audioread('../Speech_Signals/White_noise1.wav'); 
resample_noise = resample(y_noise,fs_RIR,Fs_noise);
filt_noiseL = fftfilt(RIR_noise(:,1),resample_noise);
filt_noiseL = filt_noiseL(1:sigLenSec*fs_RIR);
filt_noiseR = fftfilt(RIR_noise(:,2),resample_noise);
filt_noiseR = filt_noiseR(1:sigLenSec*fs_RIR);
filt_noise = [filt_noiseL filt_noiseR];
% Plot the noisy signal
figure(2); clf;
plot(1:size(filt_noise,1),filt_noise)
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
    disp('Final SNR')
    disp(SNR_left)
    disp('Signal scaling')
    disp(0.9999^count)
    sig_scale = 0.9999^count;
    binaural_sig_scaled = [convL*sig_scale, convR*sig_scale];
end
%% MFxLMS

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

%%
% Calculate the noise suppression
NSPL = 10*log10(mean(e_L(1,:).^2)/mean(d(:,1).^2))
NSPR = 10*log10(mean(e_R(1,:).^2)/mean(d(:,2).^2))

%% Plotting Fx-NLMS outputs
figure
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
figure
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



