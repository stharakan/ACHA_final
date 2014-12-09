% A_t.m
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

function b = A_t(x, OMEGA, P, c, flag)



N = length(x);
k = length(OMEGA);

b=toeplitzmult(c,x);
switch flag
    case 'RPMS'
        b = rpms(b,P,OMEGA);
    case 'SUB'
        b = b(OMEGA)*sqrt(N);
end