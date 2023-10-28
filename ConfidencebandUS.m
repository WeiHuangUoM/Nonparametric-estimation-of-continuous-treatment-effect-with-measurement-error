%Estimate Confidence intervals using undersmoothing

%Sub-routines: (1) CV_fdecCond.m
%              (2) hSIMEXUknown.m
%              (3) underSmooth.m
N = length(Y);

h21= CV_fdecCond(N,t0,Si,[X,Y],errortype,b,hPI);
h22 = CV_fdecCond(N,t0,Si,X,errortype,b,hPI);

d = size(X,2);
parfor dim = 1:d
    dim
[bwXT(dim),rhoXT(dim),~,~]=hSIMEXUknown(Si,X(:,dim),errortype,b,20);
end
[bwYT,rhoYT,~,~]=hSIMEXUknown(Si,Y,errortype,b,20);

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
[bwYtX,rhoYtX,~,~]=hSIMEXUknown(Si,Y-X*betahat,errortype,b,20);
ghatt = NWDecUknown(t0,Si,Y-X*betahat,errortype,b,bwYtX,rhoYtX);
EYtX = X*betahat + ghatt;

[bw_us,biashat,sdIF,muhat] = underSmooth(N,t0,X,Si,Y,hPI,K,...
    bwD,ridgeD,EYtX,errortype,varU,b,h21,h22,0);


CIlbest = muhat- 1.96*sdIF;
CIubest = muhat+1.96*sdIF;

%filename = sprintf('%s_%s.mat',filename,date);
%save(filename)