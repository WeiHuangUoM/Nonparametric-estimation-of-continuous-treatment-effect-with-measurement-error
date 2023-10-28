function sdIF = VarianceIF(N,t0,Si,Y,muhat,h,bw,weight,EYtX,errortype,b,rho,hPI)
% Wei Huang and Zheng Zhang (2022).
% Nonparametric Estimation of the Continuous Treatment
% Effect with Measurement Error
% Estimate the asymptotic variance

%Sub-routines: (1) get_weightCon.m
%              (2) kerU_d.m
%              (3) NWDecUknown.m
%              (4) fdecUknown.m


LUh0 = kerU_d(t0,Si,errortype,b,h,0)*N;
LUh = kerU_d(t0,Si,errortype,b,bw,0)*N;

%psi
piLh0 = weight.*LUh0;
psit = (muhat.*(LUh0-mean(LUh0,1))-(EYtX.*piLh0 - mean(EYtX.*piLh0,1)));

%phi
piYL = weight.*Y.*LUh;
phit =(piYL - mean(piYL,1) - muhat.*(LUh - mean(LUh,1)));

%variance
fT = fdecUknown(t0,Si,hPI,errortype,b);
fT = max(fT,rho);
eta = phit+psit;
eta_stand = eta - mean(eta,1);
v_eta = (eta_stand'*eta_stand)/(N-1);

v = v_eta./(N*fT.^2);
vv = diag(v);

sdIF = sqrt(vv);
sdIF = sdIF';
end