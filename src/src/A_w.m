% A_w.m
%
% Takes "scrambled Hadamard" measurements.
%
% Usage: b = A_w(x, OMEGA, P)
%
% x - N vector
%
% b - K vector
%
% OMEGA - K vector denoting which WH coefficients to use
%        
% P - Permutation to apply to the input vector.  Fourier coeffs of
%     x(P) are calculated.
%     Default = 1:N (no scrambling).
%

function b = A_w(x, OMEGA, P, flag, R)

N = length(x);
if (nargin < 3), P = 1:N;  end

if(any(abs(R) > 1))
    x = x(R);
else
    x = x.*R;
end

fx = fwht(x);
switch flag
    case 'RPMS'
        b = rpms(fx,P,OMEGA);
    case 'SUB'
        b = fx(OMEGA)*sqrt(N);
end
end
