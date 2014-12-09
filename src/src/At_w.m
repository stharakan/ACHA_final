% At_w.m
%
% Adjoint for "scrambled Fourier" measurements.
%
% Usage: x = At_w(b, N, OMEGA, P)
%
% b - K vector = [real part; imag part]
%
% N - length of output x
%
% OMEGA - K vector denoting which Fourier coefficients to use
%
% P - Permutation to apply to the input vector.  Fourier coeffs of
%     x(P) are embedded.
%     Default = 1:N (no scrambling).
%


function x = At_w(b, N, OMEGA, P, flag, R)

if (nargin < 4),  P = 1:N;  end

K = length(b);
fx = zeros(N,1);
switch flag
    case 'RPMS'
        fx = rpmst(b,OMEGA);
    case 'SUB'
       fx(OMEGA) = b;
end
%x = zeros(N,1);
if(any(abs(R) > 1))
    x(R) = ifwht(fx);
else
    x = ifwht(fx).*R;
end
x =x';
end
