%Author: Osama Ullah Khan,
%        Phd Student, University of Michigan-Ann Arbor.
%        Email: oukhan@umich.edu
%        Version: 1.0
%
%This code demonstrate compressive sensing example. In this
%example the signal is sparse in time domain and random samples
%are taken in frequency domain.

close all;
clear all;

sigma = 0.01;
variance = sigma.^2;
lambda = 2;

%number of samples per period
s=4;

%RF frequency
f=4e9;

%pulse repetition frequency
prf=1/30e-9;

%sampling frequency
fs=s*f;

%Total Simulation time
T=30e-9;

%t=0:1/fs:T;
t=linspace(0,T,512);

%generating pulse train
%x=pulstran(t,15e-9,'gauspuls',f,0.5);
x=pulstran(t,15e-9,'gauspuls',f,0.5);
x = x - mean(x);
size(x)
%x = ifwht(x);
size(x)

%length of the signal
N=length(x);

%Number of random observations to take
K=90;

figure;
subplot(2,1,1);
plot(t,x)
grid on;
xlabel('Time');
ylabel('Amplitude');
title(sprintf('Original Signal, UWB Pulse RF freq=%g GHz',f/1e9));

%taking Discrete time Fourier Transform of the signal
xf=fft(x);

xfmag=10*log10(abs(xf));

subplot(2,1,2);
plot(abs(xf))
grid on;
xlabel('Samples');
ylabel('Amplitude');
title('Discrete Fourier Transform of UWB pulse');


%creating dft matrix
B=hadamard(N);
Binv=inv(B);

%Selecting random rows of the DFT matrix
q=randperm(N);

%creating measurement matrix
A=B(q(1:K),:);

%taking random frequency measurements
y=(A*x');

diff = peak2peak(y);
% add noise scaled by the vector
%y = y+diff*normrnd(0,sigma,length(y),1);

% Calculating Initial guess
x0=A'*y;

%Running the recovery Algorithm
tic
xp=l1eq_pd(x0,A,[],y,1e-5);
toc




%creating dft matrix
B=dftmtx(N)*hadamard(N);
Binv=inv(B);

%Selecting random rows of the DFT matrix
q=randperm(N);

%creating measurement matrix
A=B(q(1:K),:);

%taking random frequency measurements
y=(A*x');

diff = peak2peak(y);
% add noise scaled by the vector
%y = y+diff*normrnd(0,sigma,length(y),1);

% Calculating Initial guess
x0=A'*y;

%Running the recovery Algorithm
tic
xp2=l1eq_pd(x0,A,[],y,1e-5);
toc


%creating dft matrix
B=eye(N);
Binv=inv(B);

%Selecting random rows of the DFT matrix
q=randperm(N);

%creating measurement matrix
A=B(q(1:K),:);

%taking random frequency measurements
y=(A*x');

diff = peak2peak(y);
% add noise scaled by the vector
%y = y+diff*normrnd(0,sigma,length(y),1);

% Calculating Initial guess
x0=A'*y;

%Running the recovery Algorithm
tic
xp3=l1eq_pd(x0,A,[],y,1e-5);
toc



figure;
subplot(2,2,1)
plot(t,ifwht(x))
grid on;
xlabel('Time');
ylabel('Amplitude');
title(sprintf('Original Signal, UWB Pulse RF freq=%g GHz',f/1e9));


subplot(2,2,2)
plot(t,real(ifwht(xp)),'r')
grid on;
xlabel('Time');
ylabel('Amplitude');
title(sprintf('Recovered Signal, identity sampling'));

subplot(2,2,3)
plot(t,real(ifwht(xp2)),'r')
grid on;
xlabel('Time');
ylabel('Amplitude');
title(sprintf('Recovered Signal, Fourier sampling'));

subplot(2,2,4)
plot(t,real(ifwht(xp3)),'r')
grid on;
xlabel('Time');
ylabel('Amplitude');
title(sprintf('Recovered Signal, Hadamard sampling'));


