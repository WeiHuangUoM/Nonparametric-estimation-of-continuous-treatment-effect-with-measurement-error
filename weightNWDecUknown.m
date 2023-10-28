function [y,Z, KUt0]=weightNWDecUknown(t0,Si,Y,weight,errortype,sigU,bw,rhogrid)

%estimate mu

% Inputs:
%t0: n-vector of t-values where to mu(t)
%Si: vector of contaminated data Si_1,...,Si_N
%Y: vector of data Y_1,...,Y_N
%weight: N x n matrix of pi_hat(t,X)
%bw: bandwidth h
%
%errortype: 'Lap' for Laplace errors and 'norm' for normal errors. For other error distributions, simply redefine phiU below 
%sigU: parameter of Laplace or normal errors used only to define phiU.
%
%rhogrid: (candidate for) the ridge parameter, can be a scaler or a vector. 

% Outputs:
% y: estimated mu
% Z: N x 1 response variable pi_hat*Y
% KUt0: N x n matrix of deconvolution kernel at t0 -  LU((t0-S)/h)/(N*h)
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%								WARNINGS:
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%
% The range of t-values -1 and 1 correspond to the support of phiK. 
% If you change phiK and take a kernel for which phiK is not supported on [-1,1] you have to change -1 and 1 accordingly.
%
% The phiK here must match the phiK used to compute the bandwidth (SIMEX or other).
%
% The estimator can also be computed using the Fast Fourier Transform, which is faster, but more complex. 
% See Delaigle, A. and Gijbels, I. (2007). Frequent problems in calculating integrals and optimizing objective functions: a case study in density deconvolution.   Statistics and Computing,  17,  349 - 355
% However if the grid of t-values is fine enough, the estimator can simply be computed like here without having problems with oscillations.
%
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------



% --------------------------------------------------------
% Preliminary calculations and initialisation of functions
% --------------------------------------------------------


%Default values of phiU(t)=characteristic function of the errors
%If you want to consider another error type, simply replace phiU by the characteristic function of your error type

if strcmp(errortype,'Lap')==1
	phiU=@(t) 1./(1+sigU^2*t.^2);
elseif strcmp(errortype,'norm')==1
	phiU = @(t) exp(-sigU^2*t.^2/2);
end


%phiK: Fourier transform of the kernel K. You can change this if you wish to use another kernel but make sure 
%you change the range of t-values, which should correspond to the support of phiK

phiK = @(t) (1-t.^2).^3;
lt = length(t0);
t0=reshape(t0,1,lt);

%Range of t-values (must correspond to the domain of phiK)
dt = .002;
t = (-1:dt:1)';
longt=length(t);
t=reshape(t,longt,1);


%Compute the empirical characteristic function of W (times n) at t/h: \hat\phi_W(t/h)
OO=Si*t'/bw;
csO=cos(OO);
snO=sin(OO);
clear OO;


xt=t/bw*t0;
cxt=cos(xt);
sxt=sin(xt);
clear xt;

phiUth=phiU(t/bw);
phiKU=reshape(phiK(t)./phiUth,1,longt);

%compute K_U((t0-W)/h);
KUt0=(csO.*phiKU)*cxt + (snO.*phiKU)*sxt;

N = length(Y);
KUt0=KUt0*(dt/(N*bw*2*pi));
clear cSO phiKU cxt snO sxt

%Compute response variable
Z=bsxfun(@times,weight,Y);

Den=sum(KUt0,1);
Num = sum(Z.*KUt0,1);

lrho = length(rhogrid);
y = zeros(lrho,lt);

%If denomintor is too small, replace it by the ridge rho
for k = 1:lrho
    ridge = rhogrid(k);
    dd = Den;
    dd(dd<ridge) = ridge;
    y(k,:) = Num./dd;
end
end
