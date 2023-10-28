function [bw,hw,ridge,idr,idx,minISE] = Opt_Calibw(t0,mu,Si,Y,X,errortype,b,K)
N = length(Y);

if strcmp(errortype,'Lap')==1
    varU = 2*b^2;
elseif strcmp(errortype,'norm')==1
    varU = b^2;
end

longbw = 10;
hPIfX=PI_deconvUknownth4(Si,errortype,varU,b);
ha=hPIfX/1.5;
hb=2*hPIfX;
bwg = ha:(hb-ha)/(longbw-1):hb;

%rho grid
hS=1.06*sqrt(var(Si))*N^(-1/5);
ab=quantile(Si,[0.05,0.95]);
xout=outerop(ab,Si,'-');
fWEF=normpdf(xout,0,hS)*ones(N,1)/N;
rhogrid=min(fWEF)*(0.1:0.025:5);

ISE= zeros(length(rhogrid),longbw);
%for j = 1:longbw
   weight = get_weightCon(t0,X,X,Si,errortype,hPIfX,K,b);
   for i =1:longbw
   [Ethyes,~,~] = weightNWDecUknown(t0,Si,Y,weight,errortype,b,bwg(i),rhogrid);
   ISE(:,i) = mean((Ethyes-mu).^2,2);
   end
%end

minISE = min(min(ISE));
minind=find(ISE==minISE, 1, 'last' );
[idr,idx]=ind2sub(size(ISE),minind);


bw = bwg(idx);
hw = hPIfX;
ridge = rhogrid(idr);
end
