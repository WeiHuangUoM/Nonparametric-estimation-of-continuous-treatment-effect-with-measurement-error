function biashat = BiasIF(N,t0,X,Si,Y,muhat,EYtX,bw,h0,errortype,varU,b,h21,h22)
% Wei Huang and Zheng Zhang (2022).
% Nonparametric Estimation of the Continuous Treatment
% Effect with Measurement Error
%Estimate bias using asymptotic theory

%Sub-routines: (1)PI_deconvUknownth4.m
%              (2)CV_fdecCond.m
%              (3)fdecCondUKnown.m
%              (4) kerU_d.m
%              (5) fdecUknown.m

kappa21 = 6; %second moment corresponding to phiK(t) = (1-t^2)^3.

hPI = PI_deconvUknownth4(Si,errortype,varU,b);

h1 = hPI;
%[h1,h21] = CV_fdecCond(N,t0,Si,[X,Y],errortype,b,hPI);
fTYX2hat = fdecCondUKnown(N,t0,Si,[X,Y],h1,h21,errortype,b,2);

%[h1,h22] = CV_fdecCond(N,t0,Si,X,errortype,b,hPI);
fTXhat = fdecCondUKnown(N,t0,Si,X,h1,h22,errortype,b,0);

Phihati = Y.*fTYX2hat./fTXhat;
Phihat = mean(Phihati,1);

fT2i = kerU_d(t0,Si,errortype,b,hPI,2)*N;
fThat = fdecUknown(t0,Si,hPI,errortype,b);

mufTi = fT2i.*muhat./fThat;
mufT = mean(mufTi,1);

biashath = kappa21/2*(Phihat - mufT)*bw^2;

fTX2hat = fdecCondUKnown(N,t0,Si,X,h1,h22,errortype,b,2);
Phitildei = EYtX.*fTX2hat./fTXhat;
Phitilde = mean(Phitildei,1);

biashath0 = kappa21/2*(mufT-Phitilde)*h0^2;

biashat = biashath + biashath0;
end

