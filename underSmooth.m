function [bw_us,biashat,sdIF,muhat] = underSmooth(N,t00,X,Si,Y,hPI,K,bwD,ridgeD,EYtX,errortype,varU,b,h21,h22,reridge)

%Choose the undersmooth bandwidth for confidence intervals

%reridge = 1 if we select new ridge parameter for each undersmooth bandwidth

%Sub-routines: (1) get_weightCon.m
%              (2) weightNWDecUknown.m
%              (3) BiasIF.m
%              (4) VarianceIF.m

weight0 = get_weightCon(t00,X,X,Si,errortype,hPI,K,b);

muhat0= weightNWDecUknown(t00,Si,Y,weight0,errortype,b,bwD,ridgeD);
biashat = BiasIF(N,t00,X,Si,Y,muhat0,EYtX,bwD,hPI,errortype,varU,b,h21,h22);
sdIF= VarianceIF(N,t00,Si,Y,muhat0,hPI,bwD,weight0,EYtX,errortype,b);

bwb=bwD;
hb=hPI/1.1;
while mean(biashat.^2) >= mean(sdIF.^2)/100 && bwb>bwD/1.98
    bwb = bwb/1.01;
    if reridge == 1
        ridgeD = SIMEXRidge(Si,Y,X,bwb,errortype,b,20,K);
    end
    muhat = weightNWDecUknown(t00,Si,Y,weight0,errortype,b,bwb,ridgeD);
    biashat = BiasIF(N,t00,X,Si,Y,muhat0,EYtX,bwb,hb,errortype,varU,b,h21,h22);
    sdIF= VarianceIF(N,t00,Si,Y,muhat,hb,bwb,weight0,EYtX,errortype,b);
end
    
bw_us = bwb;
end