function [bw,K1,K2,G] = CV_NaiveFull(Si,Y,X)

%CV for naive estimator
%sub-routines: get_weight_Naive.m; CaliEst_Naive.m

N=length(Si);

hS=1.06*sqrt(var(Si))*N^(-1/5);
bw0=hS/1.5:((hS*2-hS/1.5)/4):hS*2;

Kg = 2:4;
[Kg1,Kg2,bwg] = ndgrid(Kg,Kg,bw0);
Kg1 = Kg1(:);
Kg2 = Kg2(:);
bwg = bwg(:);
lK = length(Kg1);

w = ksdensity([X,Si],[X,Si]);

Kfold = 10;
Nvset = N/Kfold;
D = zeros(N,lK);
    for j = 1:lK
        K1 = Kg1(j);
        K2 = Kg2(j);
        bw = bwg(j);
        for i = 1:Kfold
            idx = (i-1)*Nvset+1:i*Nvset;
            D(idx,j) = diffi(idx,K1,K2,bw,X,Y,Si);
        end
    end
  
G=sum(D.*w,1);
index = find(G==min(G), 1 );
bw = bwg(index);
K1 = Kg1(index);
K2 = Kg2(index);

end

function D = diffi(i,K1,K2,bw,X,Y,Si)
Siv = Si(i);
Yv = Y(i);

Sit = Si;
Sit(i)=[];
Xt = X;
Xt(i)=[];
Y(i)=[];

weight = get_weight_Naive(Si,X,Xt,Sit,K1,K2);
weight = diag(weight);
wv = weight(i);
weight(i) = [];

[Ethyes,~] = CaliEst_Naive(bw,Siv,Sit,Y,weight);

D = (wv.*Yv - Ethyes').^2;

end

