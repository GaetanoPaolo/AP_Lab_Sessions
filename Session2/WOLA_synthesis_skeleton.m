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


L = size(X,2);

% ## Perform IFFT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section of code to complete (1 - 3 lines) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% ## Apply synthesis window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section of code to complete (1 - 3 lines) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% ## Obtain re-synthesised signals

for m = 0:L-1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Section of code to complete (1 - 3 lines) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

end
