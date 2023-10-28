function [bw,bwStar,bw1,bw2,ridge] = SIMEXRegbw(Si,Y,X,errortype,b,Bm,K)

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
%w = zeros(N,1);
%w(indw) = 1;

CritStar = zeros(longbw,length(rhogrid));
for j=1:Bm
    
    %j
    
    Wstar = Sistar(:,j);
   
    h = PI_deconvUknownth4(Wstar,errortype,varU,b)*2;
    
    try
    weight = get_weightCon(Si(indw),X,X,Wstar,errortype,h,K,b);
    catch ME
        disp(ME.message)
        h = h*1.5;
        weight = get_weightCon(Si(indw),X,X,Wstar,errortype,h,K,b);
    end
    p = diag(weight(indw,:));
    %max(p)
    
    ycv=arrayfun(@(i)weightNWDecRidgeUknown(Si(indw),Wstar,Y,weight,errortype,b,bwg(i),rhogrid),1:longbw,'UniformOutput',false);
    CV = cellfun(@(x)sum((p.*Y(indw) - x).^2,1),ycv,'UniformOutput',false);
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
    
    %j
    
    Wstar = Sidstar(:,j);
    W = Sistar(:,j);
    %weight for computing CV
    q1t=quantile(W,0.05);
    q2t = quantile(W,0.95);
    indw = W<=q2t & W>=q1t;
    clear q1t q2t

    h = PI_deconvUknownth4(Wstar,errortype,varU,b)*2;
    
    try
    weight = get_weightCon(W(indw),X,X,Wstar,errortype,h,K,b);
    catch ME
        disp(ME.message)
        h = h*1.5;
        weight = get_weightCon(W(indw),X,X,Wstar,errortype,h,K,b);
    end
    p = diag(weight(indw,:));
    %max(p)
    
    ycv=arrayfun(@(i)weightNWDecRidgeUknown(W(indw),Wstar,Y,weight,errortype,b,bwg(i),ridgej(j)),1:longbw,'UniformOutput',false);
    ycv = cell2mat(ycv);
    CV = sum((p.*Y(indw) - ycv).^2,1);
    Crit= CV';
    clear p ycv w CV
    CritStar2 = CritStar2+Crit;
    
    bw2(j) = min(bwg(Crit==min(Crit)));
end


bw = localconstant(bwStar,bw2,bw1);

ridge = mean(ridgej);
end
