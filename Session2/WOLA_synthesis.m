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
M = size(X,3);
L = size(X,2);
xs = zeros(nfft,L,M);
%x = zeros(L*
% ## Perform IFFT
for m = 1:M
    X_reconstr = cat(1,X(:,:,m),conj(flip(X(2:nfft/2,:,m),1)));
    xs(:,:,m) = ifft(X_reconstr(:,:,m),nfft);
end


% ## Apply synthesis window
for m = 1:M
%     disp(size(window))
%     disp(size(repmat(window, 1,size(xs(:,:,M),2))))
%     disp(size(xs(:,:,M)))
    
    xs(:,:,M) = xs(:,:,M).*repmat(window, 1,size(xs(:,:,M),2));
end

% ## Obtain re-synthesised signals
x = zeros((nfft/noverlap)*(L-1)+nfft,M);
for m = 1:M
    x(1:nfft,m) = xs(:,1,m);
    for l = 2:L
        x((l-1)*(nfft/noverlap)+1:(l-1)*(nfft/noverlap)+nfft,m) = x((l-1)*(nfft/noverlap)+1:(l-1)*(nfft/noverlap)+nfft,m) + xs(:,l,m);
    end
end
end

%