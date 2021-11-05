[speech, Fs] =audioread('../Speech_Signals/speech1.wav');
speech = speech(1:5*Fs);
nfft = 128;
noverlap = 2;
analwin = ones(1,(nfft/2));%+1?
synthwin = ones(1,(nfft/2));
[X,f] = WOLA_analysis(speech,Fs,analwin,nfft,noverlap);
x = WOLA_synthesis(X,synthwin,nfft,noverlap);