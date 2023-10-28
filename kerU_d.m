function y=kerU_d(t0,Si,errortype,sigU,h,l)

%Compute (the derivative of) the deconvolution kernel

%t0: a t-value where to compute the deconvoluntion kernel
%Si: vector of contaminated data Si_1,...,Si_n
%h: bandwidth
%l: order of derivative, e.g. l=1: first derivative
%
%errortype: 'Lap' for Laplace errors and 'norm' for normal errors. For other error distributions, simply redefine phiU below 
%sigU: parameter of Laplace or normal errors used only to define phiU.
%
%See
%Delaigle, A. and Hall, P. (2008). Using SIMEX for smoothing-parameter choice in errors-in-variables problems.  JASA, 103, 280-287 

%Sub-routines: outerop.m

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


%Range of t-values (must correspond to the domain of phiK)
dt = .0002;
t = (-1:dt:1)';
longt=length(t);
t=reshape(t,longt,1);

n=length(Si);

%Compute the empirical characteristic function of W (times n) at t/h: \hat\phi_W(t/h)
OO=outerop(Si,t/h,'*');
csO=cos(OO);
snO=sin(OO);
clear OO;

xt=outerop(t/h,t0,'*');
cxt=cos(xt);
sxt=sin(xt);
clear xt;

phiUth=phiU(t/h);
d_factor = (-1i*t/h).^l;
matphiKU=reshape(d_factor.*phiK(t)./phiUth,1,longt);


Num=(csO.*matphiKU)*cxt+(snO.*matphiKU)*sxt;
Num = Num*dt/(2*pi);

y=Num*(1/(n*h));
end


