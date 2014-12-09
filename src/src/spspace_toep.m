close all;
%clear all;

scale = .25;
sparsity = 1; %how much data to keep 1/sparsity*100%
lambda = 2;
snr = 5000;
type = 'RPMS'

% load some datas
%load ../data/stars512;
%x=stars512;
%load ../data/sameer512
%x=sameer512;
%load ../data/rachel379
%x=rachel;
%load ../data/oscar1078
%x=oscar;
%load ./l1magic/Data/boats.mat
%x = boats;
%load ../data/sloth500
%x = sloth500;
x = phantom('Modified Shepp-Logan',512);
x = im2double(x);
x = imresize(x,scale);
x = x/norm(x(:));
x = x - mean(x(:));
[m,n,d] = size(x);
x=x(:);

% Threshholding to increases sparsity.
s = sort(abs(x));
thresh = s(end-round(m*n*d*sparsity)+1)
x(abs(x)<thresh) = 0;
N=length(x);

% Number of observations, i.e. the number of rows in the "random" matrix
K=round(N/3);

% q is an indexing parameter that specifies which rows we use.
%q=randperm(N);

% Define the operators (matrices are super slow)
%A = @(x) dctpart(x,q(1:K));
%At = @(x) idctpart(x,q(1:K),N);


phase = zeros(2*N,1);
phase(2:N) = exp(2*pi*1i*rand(N-1,1));
phase(N+2:end) = flipud(conj(phase(2:N)));
phase(1) = 1;
phase(N+1) = 1;
c=phase;

switch type
    case 'RPMS'
        P = 1:N;
        P=P';
        OMEGA = rand(N,1);
        OMEGA(OMEGA>.5) = 1;
        OMEGA(OMEGA<=.5) = -1;
        A = @(z) A_t(z, OMEGA,K, c, 'RPMS');
        At = @(z) At_t(z, N, OMEGA,K, c, 'RPMS');
    case 'SUB'
        P = 1:N;
        P=P';
        %q=1:N;
        q = randperm(N);
        OMEGA = q(1:K)';
        %OMEGA = sort(OMEGA);
        A = @(z) A_t(z, OMEGA,K, c, 'SUB');
        At = @(z) At_t(z, N, OMEGA,K, c, 'SUB');
end

% Take dat measurement
y = A(x);
% add noise scaled by the vector
sigma = norm(y)/(2*snr*sqrt(length(y)));
e = normrnd(0,sigma,length(y),1);
y = y+e;
% First guess... maybe could also be 0?
x0 = At(y);
eps = sqrt(sigma.^2*(n+lambda*sqrt(2*n))); %use this one for gaussian noise
%eps = sigma*round(sqrt(m));
%epsilon = (2*max(abs(y))/10)/sqrt(12)*sqrt(K);


% Run the l1-magic optimizer
tic
%xp=l1eq_pd(x0,A,At,y,1e-5);
%xp = l1qc_logbarrier(x0, A, At, y, eps, 1e-1, 2, 1e-8, 600);
xp = tvqc_logbarrier(x0, A, At, y, eps, 1e-1, 5, 1e-8, 200);
toc

figure;
subplot(1,2,1)
imshow(reshape(x,m,n,d),[]);
title('Original Image');

subplot(1,2,2)
imshow(reshape(xp,m,n,d),[]);
title('Recovered Image');

norm(xp - x(:))/norm(x(:))


