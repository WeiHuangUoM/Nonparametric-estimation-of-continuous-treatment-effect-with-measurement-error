function h2 = CV_fdecCond(N,t0,Si,Z,errortype,sigU,hPI)

% Z: N*d
% Cross-validation for estimating conditional density f_T|Z(t|Z)
% Sub-routines: fdecCondUKnown.m

%Reshape
t0 = reshape(t0,1,length(t0));
Si = reshape(Si,N,1);
dt = t0(2)-t0(1);

%Znorm = sqrt(sum(Z.^2,2));
sdZ = std(Z,1);

ha2 = 1.06*min(sdZ)*N^(-1/5)/2;
hb2 = 1.06*max(sdZ)*N^(-1/5)*3;

hc2 = ha2:(hb2-ha2)/9:hb2;

hc1 = hPI;

[hc1,hc2]=ndgrid(hc1,hc2);
hc1 = hc1(:);
hc2 = hc2(:);

lb = length(hc1);

LUS = IbFourier(Si,Si,errortype,sigU,hPI);
LUS(1:1+size(LUS,1):end)=0;

Kfold = 10;
Kiter = N/Kfold;
    D = zeros(N,lb);
    for j = 1:lb
        [~,KerZ,LU] = fdecCondUKnown(N,t0,Si,Z,hc1(j),hc2(j),errortype,sigU,0);
        
        for k = 1:Kiter
            index = (k-1)*Kfold+1:k*Kfold;
            D(index,j) = diffi(index,KerZ,LU,LUS,dt);
        end
    end
  
G=mean(D,1);
index=find(G==min(G), 1 );
h2 = hc2(index);

end

function D = diffi(index,KerZ,LU,LUS,dt)

KerZ(index,:) = 0;
fcv = KerZ'*LU./((sum(KerZ,1))');
fcv = fcv(index,:);
I = LUS.^0;
I(1:1+size(I,1):end)=0;
fcvS = KerZ'*LUS*KerZ./(KerZ'*I*KerZ);
fcvS = fcvS(index,index);
fcvS = diag(fcvS);

D = sum(fcv.^2,2)*dt - 2*fcvS;

end

function y=IbFourier(t0,Si,errortype,sigU,h)

%Compute 1/(2pi)*\int exp{it(Sj-Sk)} \phi_L(t h_{PI})/|\phi_U(t)|^2dt

%t0: a t-value where to compute the deconvoluntion kernel
%Si: vector of contaminated data Si_1,...,Si_n
%h: bandwidth
%
%errortype: 'Lap' for Laplace errors and 'norm' for normal errors. For other error distributions, simply redefine phiU below 
%sigU: parameter of Laplace or normal errors used only to define phiU.
%
%See
%Meister(2006) page 66.

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
d_factor = (-1i*t/h);
matphiKU=reshape(d_factor.*phiK(t)./phiUth.^2,1,longt);


Num=(csO.*matphiKU)*cxt+(snO.*matphiKU)*sxt;
Num = Num*dt/(2*pi);

y=Num*(1/h);
end


