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

%% MFxLMS

M = 400;  % Length of secondary path (RIR)
L = 400;  % Adaptive filter length

mu = 0.5;   % Step size
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
hold off
legend('d(n)','e(n)')

figure
hold on 
plot(1:sigLenSample,filt_noise(:,2))
plot(1:sigLenSample,e_R)
hold off
legend('d(n)','e(n)')
