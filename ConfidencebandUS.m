%Estimate Confidence intervals using undersmoothing

%Sub-routines: (1) CV_fdecCond.m
%              (2) hSIMEXUknown.m
%              (3) underSmooth.m
%              (4) NWDecUknown
%              (5) EYTXLiang

j = 1;
X = XDG(:,j);
Si = SDG(:,j);
Y = YDG(:,j);
h21= CV_fdecCond(N,t0,Si,[X,Y],errortype,b,hPI(j));
h22 = CV_fdecCond(N,t0,Si,X,errortype,b,hPI(j));
[bwXT,rhoXT,~,~]=hSIMEXUknown(Si,X,errortype,b,20);
[bwYT,rhoYT,~,~]=hSIMEXUknown(Si,Y,errortype,b,20);

EXTt =  NWDecUknown(Si,Si,X,errortype,b,bwXT,rhoXT);
EYTt =  NWDecUknown(Si,Si,Y,errortype,b,bwYT,rhoYT);
Xtilde = X - EXTt';
Ytilde = Y - EYTt';
XYtilde = Xtilde'*Ytilde;
betahat = (Xtilde'*Xtilde)\XYtilde;

[bwYtX,rhoYtX,~,~]=hSIMEXUknown(Si,Y-X*betahat,errortype,b,20);

parfor j = 1:J
    j
    X = XDG(:,j);
    Si = SDG(:,j);
    Y = YDG(:,j);
    
    EYtX = EYTXLiang(t0,X,Si,Y,bwYtX,rhoYtX,bwXT,rhoXT,bwYT,rhoYT,errortype,b);
    
    [bw_us(j),h_us(j),biashat(j,:),sdIF(j,:),muhat(j,:)] = underSmooth(N,t0,X,Si,Y,hPI(j),K(j),...
        bwSIMEX(j),ridgeSIMEX(j),EYtX,errortype,varU,b,h21,h22,0);
end

CIlbest = muhat - 1.96*sdIF;
CIubest = muhat+1.96*sdIF;

inCIest = mu<=CIubest & mu>=CIlbest;

CPest = mean(inCIest,1);

 filename=sprintf('DGP%d-%s-vUXr-%.2f-N%d-CIUSh-%s.mat',DGP,errortype,vUXr,N,date);
 save(filename)
 