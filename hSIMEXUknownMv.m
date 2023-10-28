function [h,rho,bw]=hSIMEXUknownMv(Si,X,Y,errortype,b)

%Author: Wei Huang modified from Aurore Delaigle to find bandwidth for bivariate NW
%estimator computed from one contaminated covariate and one non-contaminated.
%
%

%Si: vector of contaminated data S_1,...,S_n
%X: vector of non-contaminate data X_1,...,X_n
%Y: vector of data Y_1,...,Y_n
%h: bandwidth for Si
%bw: bandwidth for X
%
%errortype: 'Lap' for Laplace errors and 'norm' for normal errors. For other error distributions, simply redefine phiU below 
%b: parameter of Laplace or normal errors used only to define phiU.
%rho: ridge parameter. 


% --------------------------------------------------------
% Preliminary calculations and initialisation of functions
% --------------------------------------------------------
%Number of SIMEX samples
BB=35;

n=length(Si);
SiS = repmat(Si,1,BB);
if strcmp(errortype,'Lap')==1
	Ustar = (rlap(b,2*BB,n))';
    varU = 2*b^2;
elseif strcmp(errortype,'norm')==1
	Ustar = normrnd(0,b,[n,2*BB]);
    varU = b^2;
end

%Define a grid where to search for the SIMEX bandwidth for Si. By default we take [h/2,2h], where h=PI bandwidth for density estimation.
%Increase the gird if too small
hPIfX=PI_deconvUknownth4(Si,errortype,varU,b);
ha=hPIfX/1.5;
hb=2*hPIfX;
gridh=ha:(hb-ha)/19:hb;
clear ha hb

%Define a grid where to search for the SIMEX bandwidth for X. By default we take [bw/2,2bw], where bw=ROT bandwidth for density estimation.
%Increase the gird if too small
hX = 1.06*std(X)*n^(-1/5);
gridbw=hX;

%Define a grid where to search for rho. 
%Recall that rho prevents the denominator of the NW estimator from being too small. In the SIMEX world, the denominator estimates the contaminated density f_W
%This is what motivates the default grid for rho used here.

%Estimator of fW(q_{0.05}) and fW(q_{0.95}) using standard (error-free) KDE and normal reference bandwidth, where q_{alpha} denotes the alpha empirical quantile of the W_i's.
hW=1.06*sqrt(var(Si))*n^(-1/5);
ab=quantile(Si,[0.05,0.95]);
xout=outerop(ab,Si,'-');
fWEF=normpdf(xout,0,hW)*ones(n,1)/n;
gridrho=min(fWEF)*(0.025:0.025:4);
clear hW ab xout fWEF


lh=length(gridh);
	
h1 = zeros(BB,1);
rho = zeros(BB,1);
%---------------------------------------------------------------------
%Step 1: find the ridge parameter using only the first level of SIMEX
%---------------------------------------------------------------------
Sistar = SiS + Ustar(:,1:BB);

for bb=1:BB

    Wstar = Sistar(:,bb);

	%For each h in the grid of h-candidates, compute the CV criterion for the data Wstar (this will automatically consider all rho candiates)
	CV = arrayfun(@(i)MvNWDecL1OCUknown(Si,Wstar,X,Y,errortype,b,gridh(i),gridbw,gridrho),1:lh,'UniformOutput',false);
    CV = cell2mat(CV');
	%CVrho=CVrho+cell2mat(CV');
	minCV=find(CV==min(min(CV)), 1, 'last' );
    [indh,indrho]=ind2sub(size(CV),minCV);

    rho(bb)=gridrho(indrho);
    h1(bb)=gridh(indh);
    
end

%----------------------------------------
%Step 2: Keep rho fixed and find h SIMEX 
%----------------------------------------

%CVhstar=0*gridh;
h2 = zeros(BB,1);

Sidstar = Sistar + Ustar(:,BB+1:2*BB);
for bb=1:BB

	Wstar = Sidstar(:,bb);
    W = Sistar(:,bb);
	
	%Compute CV for each h in the grid, using the ridge parameter rho found above
	CV = arrayfun(@(i)MvNWDecL1OCUknown(W,Wstar,X,Y,errortype,b,gridh(i),gridbw,rho(bb)),1:lh);
    %CV = cell2mat(CV');
	%CVhstar=CVhstar+cell2mat(CV');
    
    indh=CV==min(CV);
    h2(bb)=min(gridh(indh));
    
end



h = locallinear(h1,h2,h1);
h = mean(h);
if (h==Inf || h==-Inf || isnan(h))
    h=hPIfX;
end
rho = mean(rho);
bw = gridbw;
end