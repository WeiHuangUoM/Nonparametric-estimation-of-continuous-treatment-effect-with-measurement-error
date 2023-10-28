function EYtX = EYTXLiang(t0,X,Si,Y,bwYtX,rhoYtX,bwXT,rhoXT,bwYT,rhoYT,errortype,b)
% Wei Huang and Zheng Zhang (2022).
% Nonparametric Estimation of the Continuous Treatment
% Effect with Measurement Error
% Estimate the asymptotic variance

%Sub-routines: (1) get_weightCon.m
%              (2) kerU_d.m
%              (3) NWDecUknown.m
%              (4) fdecUknown.m


%Estimate E(Yi|Ti=t,Xi=Xi) using Liang (2000)
d = size(X,2);
if d>1
    parfor dim = 1:d
        EXTt(dim,:) =  NWDecUknown(Si,Si,X(:,dim),errortype,b,bwXT(dim),rhoXT(dim));
    end
else
EXTt =  NWDecUknown(Si,Si,X,errortype,b,bwXT,rhoXT);
end
EYTt =  NWDecUknown(Si,Si,Y,errortype,b,bwYT,rhoYT);
Xtilde = X - EXTt';
Ytilde = Y - EYTt';
XYtilde = Xtilde'*Ytilde;
betahat = (Xtilde'*Xtilde)\XYtilde;
ghatt = NWDecUknown(t0,Si,Y-X*betahat,errortype,b,bwYtX,rhoYtX);

EYtX = X*betahat + ghatt;
end