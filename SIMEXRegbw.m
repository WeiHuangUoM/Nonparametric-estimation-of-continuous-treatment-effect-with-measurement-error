function [bw,bwStar,bw1,bw2,ridge] = SIMEXRegbw(Si,Y,X,K,errortype,b,Bm)

% Wei Huang and Zheng Zhang (2022).
% Nonparametric Estimation of the Continuous Treatment
% Effect with Measurement Error
% SIMEX method for choosing h to estimate mu(t)

%Inputs:
%Si: vector of contaminated data Si_1,...,Si_N
%Y: vector of data Y_1,...,Y_N
%X: N*r matrix of the data of covariables X
%K: number of sieve basis function to estimate pi_0.
%errortype: 'Lap' for Laplace errors and 'norm' for normal errors. For other error distributions, simply redefine phiU below 
%b: parameter of Laplace or normal errors used only to define phiU.
%Bm: number of bootstraps in SIMEX. Default is 35.

%Sub-routines: (1)PI_deconvUknownth4.m
%              (2) outerop.m
%              (3) get_weightCon.m
%              (4) weightNWDecRidgeUknown.m
%              (5) localconstant.m


if isempty(Bm)
    Bm=35;
end

N = length(Y);

%Generate two sets of contaminated data
SiS = repmat(Si,1,Bm);
if strcmp(errortype,'Lap')==1
	Ustar = (rlap(b,2*Bm,N))';
    varU = 2*b^2;
elseif strcmp(errortype,'norm')==1
	Ustar = normrnd(0,b,[N,2*Bm]);
    varU = b^2;
end

%bw grid
longbw = 20;
hPIfX=PI_deconvUknownth4(Si,errortype,varU,b);
ha=hPIfX/1.5;
hb=2*hPIfX;
bwg = ha:(hb-ha)/(longbw-1):hb;
clear hPIfX ha hb

%rho grid
hS=1.06*sqrt(var(Si))*N^(-1/5);
ab=quantile(Si,[0.05,0.95]);
xout=outerop(ab,Si,'-');
fWEF=normpdf(xout,0,hS)*ones(N,1)/N;
rhogrid=min(fWEF)*(0.1:0.025:5);
clear hS ab xout fWEF

Sistar = SiS + Ustar(:,1:Bm);
    
bw1 = zeros(Bm,1);
ridgej = zeros(Bm,1);

%weight for computing CV
q1t = quantile(Si,0.05);
q2t = quantile(Si,0.95);
indw = Si>q1t&Si<q2t;
w = zeros(N,1); w(indw)=1;

CritStar = zeros(longbw,length(rhogrid));
for j=1:Bm
    
    %j
    
    Wstar = Sistar(:,j);
   
    h = PI_deconvUknownth4(Wstar,errortype,varU,b);
    
    try
    weight = get_weightCon(Si,X,X,Wstar,errortype,h,K,b);
    catch ME
        disp(ME.message)
        h = h*1.5;
        weight = get_weightCon(Si,X,X,Wstar,errortype,h,K,b);
    end
    p = diag(weight);
    %max(p)
    
    ycv=arrayfun(@(i)weightNWDecRidgeUknown(Si,Wstar,Y,weight,errortype,b,bwg(i),rhogrid),1:longbw,'UniformOutput',false);
    CV = cellfun(@(x)sum((p.*Y - x).^2.*w,1),ycv,'UniformOutput',false);
    clear p ycv
    
    Crit= cell2mat(CV');
    clear CV
    
    CritStar = CritStar + Crit;
    minCV=find(Crit==min(min(Crit)) , 1, 'last' );
    [indbw,indrho]=ind2sub(size(Crit),minCV);

    %Rigdge parameter
    ridgej(j)=rhogrid(indrho);

    %h from SIMEX level 1
    bw1(j)=bwg(indbw);
end
minCV=find(CritStar==min(min(CritStar)) , 1, 'last' );
[indbw,~]=ind2sub(size(CritStar),minCV);
bwStar = bwg(indbw);
clear minCV indbw CritStar

Sidstar = Sistar + Ustar(:,Bm+1:2*Bm);
bw2 = zeros(Bm,1);
CritStar2 = zeros(longbw,1);
for j=1:Bm
    
   % j
    
    Wstar = Sidstar(:,j);
    W = Sistar(:,j);
    %weight for computing CV
    q1t=quantile(W,0.05);
    q2t = quantile(W,0.95);
    indw = W<=q2t & W>=q1t;
    w = zeros(N,1); w(indw)=1;
    
    clear q1t q2t

    h = PI_deconvUknownth4(Wstar,errortype,varU,b);
    
    try
    weight = get_weightCon(W,X,X,Wstar,errortype,h,K,b);
    catch ME
        disp(ME.message)
        h = h*1.5;
        weight = get_weightCon(W,X,X,Wstar,errortype,h,K,b);
    end
    p = diag(weight);
    %max(p)
    
    ycv=arrayfun(@(i)weightNWDecRidgeUknown(W,Wstar,Y,weight,errortype,b,bwg(i),ridgej(j)),1:longbw,'UniformOutput',false);
    ycv = cell2mat(ycv);
    CV = sum((p.*Y - ycv).^2.*w,1);
    Crit= CV';
    clear p ycv w CV
    CritStar2 = CritStar2+Crit;
    
    bw2(j) = min(bwg(Crit==min(Crit)));
end


bw = localconstant(bwStar,bw2,bw1);

ridge = mean(ridgej);
end
