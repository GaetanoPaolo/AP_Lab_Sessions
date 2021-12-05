% Session 4 - Exercise 1 - Single-channel FxNLMS ANC
% 
% Main points:
% (1) Generate the noisy microphone signal at the left ear
% (2) Implement the FxNLMS ANC.
% (3) Compute input/output SNRs and SNR improvement.
% (4) Implement the filter for left and right ears. Hear the synthesised signal

clear all;
% close all

addpath("../Session2")

% Load RIRs
load('../sim_environment/Computed_RIRs.mat')
load('../HRTF.mat') 
% Set length
J = size(RIR_sources,3);% Number of loudspeakers
Lh = 400;
sigLenSec = 5;
sigLenSample = 5*fs_RIR;
RIR = RIR_sources(1:Lh,:,:);

%%
% Plot the RIRs of the noise
figure(1); clf;
plot(1:length(RIR_noise(1:1500,1)),RIR_noise(1:1500,1))

%%

% Read in the noise source (resample if necessary)
noisefilename = {'../Speech_Signals/Babble_noise1.wav','../Speech_Signals/White_noise1.wav'};
noise_amount = size(RIR_noise,3);
[response_length,mic_amount,speaker_amount] = size(RIR_sources);
resample_noise_signals = [];
max_length=sigLenSec*fs_RIR;
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

% Plot the noisy signal
figure(2); clf;
plot(1:length(resample_noise_signals),resample_noise_signals)

%% FxLMS

M = 400;  % Length of secondary path (RIR)
L = 400;  % Adaptive filter length

mu = 0.5;   % Step size
delta = 5*10^(-5);

W = zeros(L,1); % Initialize adaptive filter
filt_noise = [zeros(L+M-1,1); resample_noise_signals];
e = zeros(sigLenSample,1);    

tic
for n = 1:sigLenSample

    % STEP 1 : Arrange the previous L + M − 1 samples of x(n) up until the
    % current sample, n, in an [M x L] Hankel matrix X_Hmat (this will be
    % use for the filtering operation)
    temp = filt_noise(n:n+(L+M-1));
    X_Hmat = hankel(temp(1:L),temp(L+1:end));

    % STEP 2: Filter the signal x(n) with the adaptive filter w(n) to
    % generate loudspeaker signals y(n). Store the samples of y(n) into the
    % vector y
    y = X_Hmat*W;
    
    % STEP 3: Compute the error signal e(n)
    e(n) = h'*y;

    % STEP 4: Pre-filter x(n) with the (estimated) RIRs to obtain the
    % filtered x(n), \bar{x}(n). Store the samples of \bar{x}(n) into
    % the vector xf
    xf = temp; %Need the correct RIRs 
    
    % STEP 5: Update the filter w(n). Store the samples of w(n) into
    % the vector w
    w(n) = w(n-1) - xf*e(n)*mu/(norm(xf,'fro')+delta)
end
toc


%%
% Calculate the noise suppression


%%
% In the existing plot of the noisy signal, superimpose the error signal to
% appreciate the noise cancellation

figure(2); hold on;
