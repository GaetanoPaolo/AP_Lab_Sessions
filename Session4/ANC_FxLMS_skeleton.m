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
load('../sim_environment/Computed_RIRs_session4.mat');
% Set length
sigLenSec =10;

%%
% Plot the RIRs of the noise
[len,~] = size(RIR_noise);
figure(1); clf;

subplot(2,1,1)
plot(1:len,RIR_noise(:,1))
subplot(2,1,2)
plot(1:len,RIR_noise(:,2))
%% Generate noisy mic signal

% Read in the noise source (resample if necessary)
[y_noise,Fs_noise] = audioread('../Speech_Signals/White_noise1.wav'); 
resample_noise = resample(y_noise,fs_RIR,Fs_noise);
filt_noise = fftfilt(RIR_noise(:,1),resample_noise);
filt_noise = filt_noise(1:sigLenSec*fs_RIR);
% Plot the noisy signal
figure(2); clf;
plot(1:size(filt_noise,1),filt_noise)


%% FxLMS

M = 70;  % Length of secondary path (RIR)
L = 600;  % Adaptive filter length

mu = 0.14;   % Step size
delta = 5*10^(-5);

W = zeros(L,1); % Initialize adaptive filter

sigLenSample = sigLenSec*fs_RIR;
x = cat(1,zeros(L+M-1,1),resample_noise(1:sigLenSample));
e = zeros(1,sigLenSample);
h = RIR_sources(1:M,1,4);
d = filt_noise;

tic
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
    
    % STEP 5: Update the filter w(n). Store the samples of w(n) into
    % the vector w

    W = W - (mu/(norm(xf)^2+delta))*xf*e(n);
end
toc
%% Plotting Fx-NLMS outputs
figure
hold on 
plot(1:sigLenSample,filt_noise)
plot(1:sigLenSample,e)
hold off
legend('d(n)','e(n)')

%%
% Calculate the noise suppression
NSP = 10*log10(mean(e.^2)/mean(d.^2))

%% 
% In the existing plot of the noisy signal, superimpose the error signal to
% appreciate the noise cancellation
%4.1.7: The error e(n) doesnt go to zero, there is too much variability in the
%noise to be cancelled

figure(2); hold on;
