[speech, Fs] =audioread('../Speech_Signals/speech1.wav');
speech = speech(1:5*Fs);
nfft = 128;
noverlap = 2;
dirac = zeros(1,(nfft/2));
dirac(1) = 1;
[X,f] = WOLA_analysis(speech,Fs,dirac,nfft,noverlap);
