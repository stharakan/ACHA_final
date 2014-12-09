% A_c.m
%
% Takes "scrambled convolution" measurements.
%
% Usage: b = A_c(x, OMEGA, P, phase)
%
% x - N vector
%
% b - K vector
%
% OMEGA - K vector denoting which conv coefficients to use
%        
% P -length of output vector
%
% flag - subsampling thing.
%

function b = A_c(x, OMEGA, P, phase, flag)



N = length(x);
k = length(OMEGA);

fx = fft(x);
fx = fx.*phase;
b = ifft(fx);
switch flag
    case 'RPMS'
        b = rpms(b,P,OMEGA);
    case 'SUB'
        b = b(OMEGA)*sqrt(N);
end