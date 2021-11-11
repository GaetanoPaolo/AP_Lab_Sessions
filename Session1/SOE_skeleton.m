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
sweetspotFlag = 1;
% Flag for OLA testing time comparison
OLA_time = 0;
% Number of loudspeakers
J = size(RIR_sources,3);

% Define the lengths of RIRs and g filters
Lh = 400; % Length of impulse responses of listening room
Lg = 2*(Lh-1)/(J-2);      % Length of filters g_j

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
if noiseFlag == 0
    % Without noise
    g = H_1\x_1;
else
    % With noise
    dev= 0.05*std(H_1(:,1));
    H_1_noise = H_1+ dev*randn(size(H_1));
    g = H_1_noise\x_1;
end

%save('g.mat','g');

% Plot estimated and real HRTFs
figure(1);
hold on
plot(1:size(x_1,1),x_1);
plot(1:size(x_1,1),H_1*g);
legend('x','H*g');
title('Approximated and real HRTFs')


% Calculate synthesis error
synth_error=norm(H_1*g-x_1);
disp('Norm of error between given and synthetized HRTFs')
disp(synth_error);

% Synthethize the binaural speech using H and g and compare it
% (psychoacoustically) to the binaural speech synthetized with x

% importing speech signal
speech1 = audioread('../Speech_Signals/speech1.wav');
speech1_res = resample(speech1,8000,44100);
speech1_cut = speech1_res(1:10*8000);

%wiith standard conv
transf_tot = H_1*g;
transf_L = transf_tot(1:length(transf_tot)/2,:);
transf_R = transf_tot(length(transf_tot)/2+1:end,:);
convL = conv(speech1_cut,transf_L);
convR = conv(speech1_cut,transf_R);
binaural_sig = [convL convR];

%with OLA
nfft = 512;
convL_x_OLA = OLA(speech1_cut,transf_L,nfft);
convR_x_OLA = OLA(speech1_cut,transf_R,nfft);
binaural_sig_OLA = [convL_x_OLA,convR_x_OLA];

%with WOLA
noverlap = 2;
analwin = sqrt(hann(nfft,'periodic'));
synthwin = sqrt(hann(nfft,'periodic')); 
channels = 5;
speeches = repmat(speech1_cut,1,channels);
[X,f] = WOLA_analysis(speeches,fs_RIR,analwin,nfft,noverlap,g);
x = WOLA_synthesis(X,synthwin,nfft,noverlap);
x = repmat(x,1,2);
x_left = x(:,1:channels);
x_right = x(:,channels+1:channels*2);
synth_yeet = [x_left(:,1) x_right(:,1)];
for j = 1:channels
    x_left(:,j) = OLA(x_left(:,j),RIR_sources(1:400,1,j),nfft);
end
for j = 1:channels
    x_right(:,j) = OLA(x_right(:,j),RIR_sources(1:400,2,j),nfft);
end
binaural_sig_WOLA = [sum(x_left,2) sum(x_right,2)];

%plotting the signals
time = min([size(binaural_sig,1),size(binaural_sig_OLA,1),size(binaural_sig_WOLA,1)]);
figure
subplot(2,1,1)
title('Different filter methods for left ear')
hold on
plot(1:time,binaural_sig_WOLA(1:time,1))
plot(1:time,binaural_sig(1:time,1))
plot(1:time,binaural_sig_OLA(1:time,1))
hold off
xlabel('timesamples')
ylabel('Signal amplitude')
legend('WOLA','time domain convolution','OLA')
subplot(2,1,2)
title('Different filter methods for right ear')
hold on
plot(1:time,binaural_sig_WOLA(1:time,2))
plot(1:time,binaural_sig(1:time,2))
plot(1:time,binaural_sig_OLA(1:time,2))
hold off
legend('WOLA','time domain convolution','OLA')
xlabel('timesamples')
ylabel('Signal amplitude')
%(Synthetizing with the original HRTF's contained in X (rhs of SOE))
xL_1 =cat(1,zeros(3,1),x_1(1:(length(transf_tot)/2)-3,1));
xR_1 = x_1((length(transf_tot)/2)+1:end,1);
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

%Norm of the error between generated binaural sigs
normerror = norm(binaural_sig-binaural_synth_x);
disp('Absolute error between given and synthetized HRTFs for generated binaural sigs')
disp(normerror);

% Checking the sweetspot for 
synth_errors =zeros(1,10);
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
            temp_mat = toeplitz(RIR_sources_new(1:Lh,1,j),zeros(Lg,1));
            HL_new = [HL_new temp_mat];
        end
        for j = 1:speaker_new
            temp_mat = toeplitz(RIR_sources_new(1:Lh,2,j),zeros(Lg,1));
            HR_new = [HR_new temp_mat];
        end
        H_new = cat(1,HL_new,HR_new);
        idx_nonzeroslines = sum(abs(H_new),2)> 0;
        H_1_new = H_new(idx_nonzeroslines,:);
        transf_tot_new = H_1_new*g;
        synth_errors(i) = norm(transf_tot_new(1:738)-x_1);
    end
    figure
    plot(1:10,synth_errors)
    title('synth error between synth and original HRTFs')
    xlabel('Horiz shift of mics [cm]')
    ylabel('Norm')
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


binaurals_sig1=binaurals_sig1_1+binaurals_sig1_2;
binaurals_sig2=binaurals_sig2_1+binaurals_sig2_2;
binaurals_sig3=binaurals_sig3_1+binaurals_sig3_2;
binaurals_sig4=binaurals_sig4_1+binaurals_sig4_2;
%soundsc(binaurals_sig4_2, f_res)






