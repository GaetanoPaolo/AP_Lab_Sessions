% Lab 2 for Digital Audio Signal Processing Lab Sessions
% Session 2: Binaural synthesis and 3D audio: OLA and WOLA frameworks
% R.Ali, G. Bernardi, J.Schott, A. Bertrand
% 2021
%
% The following is the skeleton code for the OLA function, which you need to
% complete.

function y = OLA(x,h,nfft)
%
% Overlap and add method for computing convolution between the filter, h,
% and signal x.
% INPUT:
%   x           : input time-domain signal(s) (samples x 1)
%   h           : filter (samples x 1)
%   nfft        : FFT size
%
% OUTPUT:
%   y           : convolved output signal

narginchk(2,3)

Lh = length(h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section of code to complete (1-2) lines %
% Calculate the appropriate value of Lx, the frame length of x,
% given your nfft and Lh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = nfft - Lh +1;
H = fft(h, nfft);   % DFT of h
H = H(:);           % make into a column vector (speed)
H_matrix = diag(H);
nx = length(x);     % total length of x
y = zeros(nx,1);

istart = 1;
x_frame = zeros(1,nfft);
while istart <= nx
    disp(istart);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Section of code to complete (5 - 10 lines) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Lx+istart-1 > nx 
        rem_len = nx - istart+1;
        x_frame = zeros(1,nfft);
        x_frame(1:rem_len) = x(istart:istart +rem_len-1);
        x_fft = fft(x_frame,nfft);
        disp(size(ifft(H_matrix*x_fft')));
        disp(size(y(istart:nx)));
        temp_frame = ifft(H_matrix*x_fft');
        y(istart:nx) = y(istart:nx)+temp_frame(1:rem_len);
    else
        x_frame(1:Lx) = x(istart:istart +Lx-1);
        x_fft = fft(x_frame,nfft);
        disp(istart+nfft-1);
        if istart+nfft-1>nx 
            temp_frame = ifft(H_matrix*x_fft');
            y(istart:nx) = y(istart:nx)+temp_frame(1:nx-istart+1);
        else
            y(istart:istart +nfft-1) = y(istart:istart +nfft-1)+ifft(H_matrix*x_fft');
        end
    end
    
%     disp(size(H_matrix))
%     disp(size(x_fft))
%     disp(size(ifft(H_matrix*x_fft')))
%     disp(size(y(istart:istart +nfft-1)))
    
    
    
    istart = istart +Lx;
end


end
