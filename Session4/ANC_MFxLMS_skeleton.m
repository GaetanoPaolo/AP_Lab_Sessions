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

% Load RIRs

% Set length
% sigLenSec = ;


% Run the cross cancel from session 1 to obtain the auralized speech

%% Plot the RIRs of the noise
figure(1); clf; 

%%

% Read in the noise source (resample if necessary)

% Boolean flag to decide whether to add or not thee auralized speec
% addBinSig = 1;

% Plot the noisy signals (1 for each subplot)
figure(2); clf;


%% MFxLMS

% M = 400;  % Length of secondary path (RIR)
% L = 400;  % Adaptive filter length

% mu = 0.5;   % Step size

% W = zeros(L,J); % Initialize adaptive filter

tic
for n=1:sigLenSample
    
    % STEP 1 : Arrange the previous L + M − 1 samples of x(n) up until the
    % current sample, n, in an [M x L] Hankel matrix X_Hmat (this will be
    % use for the filtering operation). Store the samples in the
    % appropriate form.

    % STEP 2: Filter the signal x(n) with the adaptive filter w(n) to
    % generate loudspeaker signals y(n). Store the samples in
    % the appropriate form.
    
    % STEP 3: Compute the error signals e_L(n) and e_R(n). Store them in
    % the appropriate form.
    
    % STEP 4: Pre-filter x(n) with the (estimated) RIRs to obtain the
    % filtered x(n), \bar{x}(n). Store the samples in the appropriate form.
    
    % STEP 5: Update the filter w(n). Store the samples in the appropriate
    % form.
    
end
toc

%%
% Calculate the noise suppression

%%
% In the existing plot of the noisy signals, superimpose the corresponding
% error signals to appreciate the noise cancellation

figure(2); 

