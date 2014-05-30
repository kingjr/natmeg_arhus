function [at,W,WT]=whwt(x,f,sigma,fs)
% [x_hwt,W,WT]=whwt(x,f,sigma,fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply windowed harmonic wavelet transform of a signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       - input:
%           x: signal 
%           f: frequency of interest
%           sigma: window size
%           fs: sampling frequency
%       - output:
%           at: "The a(t) is composed of a real and an imaginary part, and the
%               imaginary part is the quadrature component of the real part
%               in theory the WHWT can be used to estimate the instantaneous
%               phase of a signal within the expected frequency band economically. 
%               Through the WHWT, the signal can be arbitrarily interpolated 
%               along the frequency axis with a constant space along the
%               time axis, so the instantaneous phase of the signal at a 
%               frequency band can be estimated."
%           W: wavelet filter used in frequency domain
%           WT: wavelet filter in time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adapted from Li et al, Neurocomputing 2011, 3389-3403
% and Li Duan scipt:
% available http://www.mathworks.com/matlabcentral/fileexchange/25727-hwtps
% (c) JeanRemi King
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. parameters
% 1.1 signal parameters
L       = length(x);
df      = fs/L;
N       = ceil(fs/df);% number of points for the fft
% 1.2 wavelet parameters
A       = 1/sqrt(3);% window parameter by default
B       = 1/sqrt(3);
s       = pi*sqrt(1/3+(B^2-16*A*B)/(4*A^2+2*B^2)/pi^2)/sigma;
w0      = f*s;% central frequency in Hz
w       = [w0/s-1/2/s+df:df:w0/s+1/2/s];

% 2. compute wavelet filter in spectral domain
W=zeros(1,N); % initialize
for i=fix(w/df)+1 % for each frequency
    W(i)=A+B*cos(2*pi*s*((i-1)*df-w0/s));
end
% 
% 3. compute wavelet filter in time domain:
% the WHWT of the signal x(t) denoted by a(t) can be obtained by taking the 
% inverse Fourier transform of A(f).
t       =[-4+1/fs:1/fs:4]; % time
w       =ifft(W,length(t)); % 
WT      =[conj(w(4*fs:-1:1)),w(1:4*fs)];

% filter signal
% 1. apply fourier transform of signal x 
Xf      = fft(x,N); % N points FFT of the signal

% 2. the multiplication of the X(f) and the conjugate of the windowed
% harmonic wavelet Ww(f) in the frequency deomain is A(f)=X(f)Ww*(f)
Af=Xf.*conj(W); 

% 3. the WHWT of the signal x(t) denoted by a(t) can be obtained by taking 
% the inverse Fourier transform of A(f)
at=ifft(Af);