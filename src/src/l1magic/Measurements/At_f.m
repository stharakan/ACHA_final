% At_f.m
%
% Adjoint for "scrambled Fourier" measurements.
%
% Usage: x = At_f(b, N, OMEGA, P)
%
% b - K vector = [real part; imag part]
%
% N - length of output x
%
% OMEGA - K/2 vector denoting which Fourier coefficients to use
%         (the real and imag parts of each freq are kept).
%
% P - Permutation to apply to the input vector.  Fourier coeffs of
%     x(P) are embedded.
%     Default = 1:N (no scrambling).
%
% Written by: Justin Romberg, Caltech
% Created: October 2005
% Email: jrom@acm.caltech.edu
%


function x = At_f(b, N, OMEGA, P, flag, R)

if (nargin < 4),  P = 1:N;  end

K = length(b);
fx = zeros(N,1);


switch flag
    case 'RPMS'
        b = 1/sqrt(2)*b(1:K/2) + 1i/sqrt(2)*b(K/2+1:K);
        fx = rpmst(b,OMEGA(1:N/2));
        fx = [fx; -fx(1);flipud(conj(fx(2:N/2)))];
    case 'SUB'
        fx(OMEGA) = sqrt(2)*b(1:K/2) + 1i*sqrt(2)*b(K/2+1:K);
end

x = zeros(N,1);

if(any(abs(R) > 1))
    x(R) = sqrt(N)*real(ifft(fx));
else
    x = sqrt(N)*real(ifft(fx)).*R;
end
x =x;
