clear all
close all
clc

% cd 'C:\Users\anchi\Desktop\hadamard\New folder'
load sigValue_randompatten4
load randompattern4

H0 = reshape(H1,4096,4096);

clear H1;
Npix=64^2; % number of (super)pixels in Hadamard pattern
Nr=64;
Nc=64;
% H0=hadamard(Npix); %generates -1s and 1s
% H0=(H0+1)/2; % to generate 0s and 1s.
% X0=imread('test.png');
% X0=imresize(double(X0),[64,64]);
% H0=(H0+1)/2; % to generate 0s and 1s.

H0=H0';
H0 = mat2gray(H0);
% Y = xlsread('GI_64x64p_4096I_3.4.19_exp2_10vG.xlsx');
% Y = xlsread('GI_64x64p_4096I_3.4.19_exp1.xlsx');

Y = sigValue(1,:);
Y=Y(:);
% H0=sparse(H0);
% X0=X0*0;
% X0(:,1:32)=1;
% X0(:,33:end)=-1;
% Y=H0*X0(:);
H1=zeros(Nr,Nc);
% H1(1,1)=1;
H1(1,1)=-4;
H1(end,1)=1;
H1(1,2)=1;
H1(2,1)=1;
H1(1,end)=1;
% H1(1,1)=-6;
% H1(end,1)=1;
% H1(1,2)=1;
% H1(2,1)=1;
% H1(1,end)=1;
% H1(2,2)=0.5;
% H1(2,end)=0.5;
% H1(end,end)=0.5;
% H1(end,2)=0.5;
N=Nr*Nc;
h=H1(:);
c=fft(h)./fft([1;zeros(N-1,1)]);

% A=H0'*H0;
lambda=1;
X=ADMM_deconv(Y(2:end),H0(2:end,:),c,lambda);
% X2 = X(10:end,10:end);
% 
figure
imagesc(X)
title('Experimental')
figure
plot(Y,'r')
hold on
plot(H0*X(:),'b')