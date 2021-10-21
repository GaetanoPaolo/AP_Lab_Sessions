%% Loading necessary stuff

load('sim_environment/Computed_RIRs.mat')
speechfilename = {'Speech_Signals/speech1.wav','Speech_Signals/speech2.wav','Speech_Signals/speech3.wav','Speech_Signals/speech4.wav','Speech_Signals/speech5.wav','Speech_Signals/speech6.wav','Speech_Signals/speech7.wav','Speech_Signals/speech8.wav'};
noisefilename = {'Speech_Signals/Babble_noise1.wav','Speech_Signals/White_noise1.wav'};
[response_length,mic_amount,speaker_amount] = size(RIR_sources);
noise_amount = size(RIR_noise,3);
if isempty(RIR_noise)
    noise_amount = 0;
end
% Params
time_seg = 10; %filtered time window in seconds
reverb = 1;
max_length=time_seg*fs_RIR;
resample_speech_signals = zeros(max_length,1);
%% plotting RIR
%disp(size(RIR_sources));
%figure
%plot(RIR_sources(:,3,2));

%% Resample speech signals
for i = 1:speaker_amount
    if (i==1) 
        [y_speech,Fs_speech] = audioread(speechfilename{i}); 
        resample_speech = resample(y_speech,fs_RIR,Fs_speech);
        resample_speech_signals = resample_speech(1:max_length,1);
    else
        [y_speech,Fs_speech] = audioread(speechfilename{i}); 
        resample_speech = resample(y_speech,fs_RIR,Fs_speech);
        resample_speech_signals = [resample_speech_signals resample_speech(1:max_length,:)];
    end
end


%% Resample noise signals
resample_noise_signals = [];
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
%% Create micsigs
% ----- Mic Generation ----- 

% Global Params:
back_noise = zeros(max_length,mic_amount);
mic = zeros(max_length,mic_amount);
% Generate noise signals in mic
for j = 1:mic_amount
    for i=1:noise_amount
        back_noise(:,j)=...
            back_noise(:,j)+...
            fftfilt(RIR_noise(:,j,i),resample_noise_signals(:,i));
    end
end

% Generate combined signals in mic
for j=1:mic_amount
    filt_speech = zeros(max_length,1);
    for i =1:speaker_amount
        filt_speech= filt_speech + fftfilt(RIR_sources(:,j,i), resample_speech_signals(:,i));
    end
    mic(:,j)=filt_speech+back_noise(:,j);
end
save('mic.mat','mic','fs_RIR');
figure
hold on
plot(1:max_length,mic(:,1));
plot(1:max_length,mic(:,2));
hold off
legend('Mic1','Mic2');