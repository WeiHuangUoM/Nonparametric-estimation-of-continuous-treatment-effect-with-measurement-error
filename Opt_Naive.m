function [bw,idx,K1,K2,minISE] = Opt_Naive(t0,mu,Si,Y,X)
N=length(Si);

%ndgrid of K1, K2
Kg = 2:4;
[Kg1,Kg2] = ndgrid(Kg,Kg);
Kg1 = Kg1(:);
Kg2 = Kg2(:);
lK = length(Kg1);

hS=1.06*sqrt(var(Si))*N^(-1/5);
bw0=hS/4:((hS*2-hS/4)/50):hS*2;
lb = length(bw0);

ISE= zeros(lb*lK,1);
for i = 1:lK
    K1 = Kg1(i);
    K2 = Kg2(i);
    weight = get_weight_Naive(Si,X,X,Si,K1,K2);
    weight=diag(weight);
    for j = 1:lb
        Ethyes= CaliEst_Naive(bw0(j),t0,Si,Y,weight);
        ISE((i-1)*lb+j) = mean((Ethyes-mu).^2);
    end
end

idx = find(ISE==min(ISE), 1 );
idxj = mod(idx,lb);
if idxj == 0
    idxj=lb;
end
idxi = (idx-idxj)/lb +1;
bw = bw0(idxj);
K1 = Kg1(idxi);
K2 = Kg2(idxi);
minISE = min(ISE);
end
