% Lab 2 for Digital Audio Signal Processing Lab Sessions
% Session 2: Binaural synthesis and 3D audio: OLA and WOLA frameworks
% R.Ali, G. Bernardi, J.Schott, A. Bertrand
% 2021
%
% The following is the skeleton code for the synthesis stage of the WOLA method, which you need to
% complete


function x = WOLA_synthesis(X,window,nfft,noverlap)
%WOLA_synthesis inverse short-time fourier transform.
%
% INPUT:
%   X           : input matrix (bins x frames x channels)
%   window      : window function
%   nfft        : FFT size
%   noverlap    : frame overlap; default: 2 (50%)
%
% OUTPUT:
%   x           : output time signal(s)

N_half = nfft / 2 + 1;
M = size(X,3);
L = size(X,2);
xs = zeros(N_half,L,M);
x = zeros(L*
% ## Perform IFFT
for m = 1:M
    xs(:,:,M) = ifft(X(:,:,M),N_half);
end


% ## Apply synthesis window
for m = 1:M
    xs(:,:,M) = xs(:,:,M).*window;
end

% ## Obtain re-synthesised signals

for m = 1:M
    for l = 0:L-1
        
    end
end

end

%