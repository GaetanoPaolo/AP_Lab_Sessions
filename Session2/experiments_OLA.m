dirac = zeros(1,50);
dirac(1) = 1;
[speech, Fs] =audioread('../Speech_Signals/speech1.wav');
speech = speech(1:5*Fs);
y = OLA(speech,dirac,100);

