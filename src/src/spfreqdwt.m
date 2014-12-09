close all;
clear all;

scale = .5;
sparsity =1; %how much data to keep sparsity*100%
snr=4.5;
lambda = 2;

% load some datas.
load ../data/sameer512
x = sameer512;
%x = phantom('Modified Shepp-Logan',512);
x = im2double(x);
x = imresize(x,scale);
%x = x/norm(x(:));
%x = x - mean(x(:));
x_orig = x;

%x = x+ imnoise(x,'gaussian',0,variance);
[m,n] = size(x);
xd = dwtf(x,'db8');
length(xd)
mt = sqrt(length(xd));
nt=mt;

%Thresholding to increase the sparsity.
s = sort(abs(xd));
thresh = s(end-round(mt*nt*sparsity)+1);
xd(abs(xd)<thresh) = 0;
N = length(xd);


% Number of observations, i.e. the number of rows in the "random" matrix
K=round(m*n/3);

figure;
subplot(1,2,1);
imshow(reshape(x,m,n),[]);

subplot(1,2,2);
imshow(abs(reshape(xd,mt,nt)),[])

% q is an indexing parameter that specifies which rows we use.
%q=randperm(N);

% Define the operators (matrices are super slow)
%A = @(x) dctpart(x,q(1:K));
%At = @(x) idctpart(x,q(1:K),N);

P = randperm(N)';
q = randperm(N/2-1)+1;
OMEGA = q(1:K/2)';

% measurement implicit matrices
A = @(z) A_f(z, OMEGA,P);
At = @(z) At_f(z, N, OMEGA,P);


% Take dat measurement
y = A(xd);
sigma = norm(y)/(2*snr*sqrt(length(y)));
e = normrnd(0,sigma,length(y),1);
% add noise scaled by the vector
y = y+e;
% First guess... maybe could also be 0?
x0=At(y);
eps = sqrt(sigma.^2*(n+lambda*sqrt(2*n))); %use this one for gaussian noise
%eps = sigma*round(sqrt(m));

% Use the l1-magic optimizer
tic
xp=l1eq_pd(x0,A,At,y,eps);
%xp = tvqc_logbarrier(x0, A, At, y, eps, 1e-1, 2, 1e-8, 600);
toc

% Put the signal back into the spatial domain
xprec=real(idwtf(xp,mt,nt,'db8'));

figure;
subplot(1,2,1);
imshow(reshape(x,m,n),[]);

subplot(1,2,2)
imshow(reshape(xprec,m,n),[])

norm(xprec - x)/norm(x);



