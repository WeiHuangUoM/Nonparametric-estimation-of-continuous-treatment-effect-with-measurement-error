function [bw_us,h_us,biashat,sdIF,muhat] = underSmooth(N,t0,X,Si,Y,h,K,bw,rho,EYtX,errortype,varU,b,h21,h22,reridge)

%Choose the undersmooth bandwidth h for confidence intervals

%reridge = 1 if we select new ridge parameter for each undersmooth bandwidth

%Sub-routines: (1) get_weightCon.m
%              (2) weightNWDecUknown.m
%              (3) BiasIF.m
%              (4) VarianceIF.m
hPI = h;
weight = get_weightCon(t0,X,X,Si,errortype,h,K,b);

[muhat,~,~] = weightNWDecUknown(t0,Si,Y,weight,errortype,b,bw,rho);
biashat = BiasIF(N,t0,X,Si,Y,muhat,EYtX,bw,h,errortype,varU,b,h21,h22);
sdIF= VarianceIF(N,t0,Si,Y,muhat,h,bw,weight,EYtX,errortype,b,rho,hPI);

muhat0=muhat;
bwb=bw;
hb =h/1.1;

while mean(biashat.^2) >= mean(sdIF.^2)/100 && bwb>bw/1.98 
    bwb = bwb/1.01;
    
    if reridge == 1
        rho = SIMEXRidge(Si,Y,X,bwb,errortype,b,20,K);
    end
    
    [muhat,~,~] = weightNWDecUknown(t0,Si,Y,weight,errortype,b,bwb,rho);
    biashat = BiasIF(N,t0,X,Si,Y,muhat0,EYtX,bwb,hb,errortype,varU,b,h21,h22);
    sdIF= VarianceIF(N,t0,Si,Y,muhat,hb,bwb,weight,EYtX,errortype,b,rho,hPI);
end
    
bw_us = bwb;
h_us=hb;
end