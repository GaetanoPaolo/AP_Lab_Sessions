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
% Flag for OLA testing time comparison
OLA_time = 0;
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
if noiseFlag == 0
    % Without noise
    g = H_1\x_1;
else
    % With noise
    dev= 0.05*std(H_1(:,1));
    %H_1 = awgn(H_1, 10*log10(40));
    H_1_noise = H_1+ dev*randn(size(H_1));
    g = H_1_noise\x_1;
end

save('g.mat','g');

% Plot estimated and real HRTFs
figure(1);
hold on
plot(1:size(x_1,1),x_1);
plot(1:size(x_1,1),H_1*g);
legend('x','H*g');


% Calculate synthesis error
synth_error=norm(H_1*g-x_1);
disp('Absolute error between given and synthetized HRTFs')
disp(synth_error);
% importing speech signal
speech1 = audioread('../Speech_Signals/speech1.wav');
speech1_res = resample(speech1,8000,44100);
speech1_cut = speech1_res(1:10*8000);
transf_tot = H_1*g;
transf_L = transf_tot(1:369,:);
transf_R = transf_tot(370:end,:);
convL = conv(speech1_cut,transf_L);
convR = conv(speech1_cut,transf_R);
binaural_sig = [convL convR];

% Synthethize the binaural speech using H and g and compare it
% (psychoacoustically) to the binaural speech synthetized with x

%Synthetizing with x
xL_1 =cat(1,zeros(3,1),x_1(1:369-3,1));
xR_1 = x_1(370:end,1);
tic
convL_x = conv(speech1_cut,xL_1);
convR_x = conv(speech1_cut,xR_1);
toc
if OLA_time
    tic
    convL_x_OLA = OLA(speech1_cut,xL_1,16384);
    convR_x_OLA = OLA(speech1_cut,xR_1,16384);
    toc
    binaural_OLA = [convL_x_OLA,convR_x_OLA];
end

    
binaural_synth_x = [convL_x convR_x];

%Syntherror
syntherror = norm(binaural_sig-binaural_synth_x);
disp('Absolute error between given and synthetized HRTFs for generated binaural sigs')
disp(syntherror);

% 1.4.10
% As the speakers are placed symetrically around the microphones, we assume
% a symmetric effect of movement and consider the right (= upward in GUI)
% movement of the speaker by steps of 10cm.
% The testing consits of computing the normalized error (like the syntherror) for each distance
% shift with the original centered microphone setup

norm_errors =zeros(1,10);
if sweetspotFlag == 1
    for i = 1:10
        %either bring the mics horizontally further away from the sources
        %m_pos(:,1) = m_pos(:,1) - ones(2,1)*0.01;
        % or shift them vertically
        m_pos(:,2) = m_pos(:,2) + ones(2,1)*0.01;
        %increasing the horizontal (=first) coordinates shortens the RIRs,
        %which brings the algorithm to failure
        [RIR_sources_new,RIR_noise_new]=create_rirs(m_pos,s_pos,v_pos,room_dim,rev_time,fs_RIR,4000);
        [~,mic,speaker_new] = size(RIR_sources_new);
        HL_new =[];
        HR_new = [];
        for j = 1:speaker_new
            temp_mat = toeplitz(RIR_sources_new(1:400,1,j),zeros(266,1));
            HL_new = [HL_new temp_mat];
        end
        for j = 1:speaker_new
            temp_mat = toeplitz(RIR_sources_new(1:400,2,j),zeros(266,1));
            HR_new = [HR_new temp_mat];
        end
        H_new = cat(1,HL_new,HR_new);
        idx_nonzeroslines = sum(abs(H_new),2)> 0;
        H_1_new = H_new(idx_nonzeroslines,:);
        transf_tot_new = H_1_new*g;
        transf_L_new = transf_tot_new(1:369,:);
        transf_R_new = transf_tot_new(370:738,:);
        convL_new = conv(speech1_cut,transf_L_new);
        convR_new = conv(speech1_cut,transf_R_new);
        binaural_sig_new = [convL_new convR_new];
        norm_errors(i) = norm(binaural_sig-binaural_sig_new);
    end
    figure
    plot(1:10,norm_errors)
    title('Normalized error between binaural signals')
    xlabel('cm')
    ylabel('abs val')
end
%% 1.4.8
% Speech1
speaker_amount = 5;
[y_speech1, f_speech1] = audioread('../Speech_Signals/speech1.wav');
f_res = 8000;
y_res = resample(y_speech1,f_res,f_speech1);
x = y_res(1:f_res*10);
binaurals_sig1_1 = [x x];
binaurals_sig2_1 = [x 0.5*x];
binaurals_sig3_1 = [x cat(1,zeros(3,1),x(1:end-3,1))];
x_left_g = fftfilt(g,x);
x_right_g = fftfilt(g,x);
filt_speech_left = zeros(80000,1);
for i =1:speaker_amount
    filt_speech_left= filt_speech_left + fftfilt(RIR_sources(:,1,i), x_left_g);
end

x_left_gh = filt_speech_left;
%filtering RIR right

filt_speech_right = zeros(80000,1);
for i =1:speaker_amount
    filt_speech_right= filt_speech_right + fftfilt(RIR_sources(:,2,i), x_right_g);
end

x_right_gh = filt_speech_right;
binaurals_sig4_1 = [x_left_gh x_right_gh];
% Speech2
[y_speech2, f_speech2] = audioread('../Speech_Signals/speech2.wav');
f_res = 8000;
y_res = resample(y_speech2,f_res,f_speech2);
x = y_res(1:f_res*10);
binaurals_sig1_2 = [x x];
binaurals_sig2_2 = [0.5*x x];
binaurals_sig3_2 = [cat(1,zeros(3,1),x(1:end-3,1)) x];
x_left_g = fftfilt(g,x);
x_right_g = fftfilt(g,x);
%filtering RIR left
filt_speech_left = zeros(80000,1);
for i =1:speaker_amount
    filt_speech_left= filt_speech_left + fftfilt(RIR_sources(:,1,i), x_left_g);
end

x_left_gh = filt_speech_left;
%filtering RIR right

filt_speech_right = zeros(80000,1);
for i =1:speaker_amount
    filt_speech_right= filt_speech_right + fftfilt(RIR_sources(:,2,i), x_right_g);
end

x_right_gh = filt_speech_right;
binaurals_sig4_2 = [x_left_gh x_right_gh];
% Copypasta
binaurals_sig1=binaurals_sig1_1+binaurals_sig1_2;
binaurals_sig2=binaurals_sig2_1+binaurals_sig2_2;
binaurals_sig3=binaurals_sig3_1+binaurals_sig3_2;
binaurals_sig4=binaurals_sig4_1+binaurals_sig4_2;
%soundsc(binaurals_sig4_2, f_res)

%% 1.4.9 
% No, due to the physical distance the soundwaves will always first arrive
% at the left ear, no matter their output





