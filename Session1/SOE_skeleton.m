% Lab 1 for Digital Audio Signal Processing Lab Sessions
% Exercise 1-4: 3D audio
% 
% In this lab, we derive a set of filters g that can be used, along with
% the measured RIRs H, to produce the proper psychocoustic rendition of 3D
% audio
%
%

clear;
% close all

% Load ATFs
load('../sim_environment/Computed_RIRs.mat')

% Load measured HRTFs
load('../HRTF.mat') 

% Define the signal length
siglength = 10;

% Load the speech signal and resample
source_filename{1} = 'speech1.wav';

% Noise flag for the noise perturbing SOE
noiseFlag = 0;
% Noise flag for sweetspot tests
sweetspotFlag = 0;

% Number of loudspeakers
J = size(RIR_sources,3);

% Define the lengths of RIRs and g filters
Lh = 400; % Length of impulse responses of listening room
% Lg = ;      % Length of filters g_j

% Truncate impulse response to reduce computational complexity

% Calculate delay for SOE

% Define the Toeplitz matrices for left and right ear (left side of SOE)
[~,mic,speaker] = size(RIR_sources);
HL =[];
HR = [];

for j = 1:speaker
    temp_mat = toeplitz(RIR_sources(1:400,1,j),zeros(266,1));
    HL = [HL temp_mat];
end
for j = 1:speaker
    temp_mat = toeplitz(RIR_sources(1:400,2,j),zeros(266,1));
    HR = [HR temp_mat];
end
        

    

% Define the HRTFs for left and right ear (right side of SOE) from the
% loaded HRTF
xL_undelayed = HRTF(:,1); % Left ear
xR_undelayed = HRTF(:,2); % Right ear
Delta=ceil(sqrt(room_dim(1)^2+room_dim(2)^2)*fs_RIR/340);
xL_delayed = cat(1,zeros(Delta,1),xL_undelayed(Delta+1:400,1));
xR_delayed = cat(1,zeros(Delta,1),xR_undelayed(Delta+1:400,1));
% Construct H (from HL and HR) and x (from xL and xR) and remove all-zero rows in H, and the corresponding elements in x
H = cat(1,HL,HR);
x = cat(1,xL_delayed,xR_delayed);
[rows,cols]=size(H);
zerows = zeros(1,rows);

idx_nonzeroslines = sum(abs(H),2)> 0;
H_1 = H(idx_nonzeroslines,:);
x_1 = x(idx_nonzeroslines,:);

% Solve the SOE
if ~noiseFlag
    % Without noise
    g = H_1\x_1;
else
    % With noise
    dev= std(H_1(:,1)); 
    g = H_1\x_1;
end

% Plot estimated and real HRTFs
figure(1);
hold on
plot(1:size(x_1,1),x_1);
plot(1:size(x_1,1),H_1*g);
legend('x','H*g');


% Calculate synthesis error
synth_error=norm(H*g-x);
disp(synth_error);
% importing speech signal
speech1 = audioread('../Speech_Signals/speech1.wav');
speech1_res = resample(speech1,8000,44100);
speech1_cut = speech1_res(1:10*8000);
transf_tot = H_1*g;
transf_L = transf_tot(1:365,:);
transf_R = transf_tot(366:end,:);
convL = conv(speech1_cut,transf_L);
convR = conv(speech1_cut,transf_R);
binaural_sig = [convL convR];

% Synthethize the binaural speech using H and g and compare it
% (psychoacoustically) to the binaural speech synthetized with x

%Synthetizing with x
xL_1 =cat(1,zeros(3,1),x_1(1:365-3,1));
xR_1 = x_1(366:end,1);
convL_x = conv(speech1_cut,xL_1);
convR_x = conv(speech1_cut,xR_1);
binaural_synth_x = [convL_x convR_x];

%Syntherror
syntherror = norm(binaural_sig-binaural_synth_x);

