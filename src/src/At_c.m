% At_c.m
%
% Adjoint for "scrambled Fourier" measurements.
%
% Usage: x = At_c(b, N, OMEGA, P)
%
% b - K vector = [real part; imag part]
%
% N - length of output x
%
% OMEGA - K vector denoting which conv coefficients to use
%
% P - length


function x = At_c(b, N, OMEGA, P, phase, flag)

fx = zeros(N,1);
switch flag
    case 'RPMS'
       fx = rpmst(b,OMEGA);
    case 'SUB'
       fx(OMEGA) = b;
end
x = zeros(N,1);
x = fft(fx);
x = x.*conj(phase);
x = ifft(x)*sqrt(N);
