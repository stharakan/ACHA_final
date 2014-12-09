close all;
clear all;

scale = 1;
sparsity = 1; %how much data to keep sparsity*100%
snr = 5000;
lambda =  2;

% load some datas.
%load ../data/sameer512
%x = sameer512;
%x = phantom('Modified Shepp-Logan',128);
idx = randperm(512*512);
idx = idx(1:round(512*512*.001));
x = zeros(512*512,1);
x(idx) = 100*rand(length(idx),1);
x = reshape(x,512,512);
x = im2double(x);
x = imresize(x,scale);
x = x/norm(x(:));
x = x - mean(x(:));
x_orig = x;

[m,n] = size(x);
%xd = dct(x(:));
xd =x(:);
%nnz(xd)
%figure()
%semilogy(sort(abs(xd),'descend'))
%xd = fft2(x);
%xd=xd(:);
%xd = fft(x(:));
length(xd)
mt = sqrt(length(xd));
nt=mt;

%Thresholding to increase the sparsity.
s = sort(abs(xd));
thresh = s(end-round(mt*nt*sparsity)+1);
xd(abs(xd)<thresh) = 0;
%x = idct(xd);
%x = idwtf(xd);
N = length(xd);

% Number of observations, i.e. the number of rows in the "random" matrix
K=round(m*n/3);

figure;
subplot(1,2,1);
imshow(reshape(x,m,n),[]);

subplot(1,2,2);
imshow(log(abs(reshape(xd,mt,nt))),[])

% q is an indexing parameter that specifies which rows we use.
%q=randperm(N);

% Define the operators (matrices are super slow)
%A = @(x) dctpart(x,q(1:K));
%At = @(x) idctpart(x,q(1:K),N);

%P = randperm(N)';
%q = randperm(N/2-1)+1;
%OMEGA = q(1:K/2)';
q = randperm(N);
P = 1:N;
P=P';
%q = 1:N;
OMEGA = q(1:K)';
%OMEGA = sort(OMEGA);
phase = zeros(m*n,1);
phase(2:m*n/2) = exp(2*pi*1i*rand(m*n/2-1,1));
phase(m*n/2+2:end) = flipud(conj(phase(2:m*n/2)));
phase(1) = 1;
phase(m*n/2+1) = 1;

% measurement implicit matrices
%A = @(z) A_f(z, OMEGA,P);
%At = @(z) At_f(z, N, OMEGA,P);
%A = @(z) A_w(z, OMEGA,P);
%At = @(z) At_w(z, N, OMEGA,P);
A = @(z) A_c(z, OMEGA,P,phase);
At = @(z) At_c(z, N, OMEGA,P,phase);

%T = eye(N);
%for i=1:100
%    %At(A(T(:,i)))
%    norm(At(A(T(:,i)))-T(:,i)*N)
%end
%A(T(:,1))



% Take dat measurement

y = A(xd);
% Generate the noise
sigma = norm(y)/(2*snr*sqrt(length(y)));
e = normrnd(0,sigma,length(y),1);
%y = y + e;
% First guess... maybe could also be 0?
x0=At(y);

eps = sqrt(sigma.^2*(n+lambda*sqrt(2*n))); %use this one for gaussian noise
%eps = sigma*round(sqrt(m));

% Use the l1-magic optimizer
tic
%xp=l1eq_pd(x0,A,At,y,1e-5);
xp = l1qc_logbarrier(x0, A, At, y, eps, 1e-1, 2, 1e-8, 600);
%xp = tvqc_logbarrier(x0, A, At, y, eps, 1e-1, 2, 1e-8, 600);
toc

% Put the signal back into the spatial domain
%xprec=real(idct(xp));
xprec = xp;
xprec = xprec/norm(xprec);
xprec = xprec - mean(xprec);

%xprec=real(ifft(xp));
%xprec=real(ifft2(reshape(xp,m,n)));

figure;
subplot(1,2,1);
imshow(reshape(x,m,n),[]);

subplot(1,2,2)
imshow(reshape(xprec,m,n),[])

norm(xprec - x(:))/norm(x(:))



