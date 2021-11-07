[speech, Fs] =audioread('../Speech_Signals/speech1.wav');
speech = speech(1:5*Fs);
nfft = 128;
noverlap = 2;
analwin = hann(nfft);%+1?
%checking 
synthwin = ones(1,nfft); 


[X,f] = WOLA_analysis(speech,Fs,analwin,nfft,noverlap);
x = WOLA_synthesis(X,synthwin,nfft,noverlap);
x_res = reshape(x',[1,size(x,1)*size(x,2)]);