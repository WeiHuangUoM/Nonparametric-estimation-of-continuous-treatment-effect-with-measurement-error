function [bw,K1,K2,G] = CV_Naive(Si,Y,X)
N=length(Si);

%ndgrid of K1, K2
Kg = 2:2;
[Kg1,Kg2] = ndgrid(Kg,Kg);
Kg1 = Kg1(:);
Kg2 = Kg2(:);
lK = length(Kg1);

%grid for bw
hS=1.06*sqrt(var(Si))*N^(-1/5);
bw0=hS/4:((hS*2-hS/4)/50):hS*2;
lb = length(bw0);

%weight
q1t=quantile(Si,0.05);
q2t = quantile(Si,0.95);
q1x = quantile(X,0.05);
q2x=quantile(X,0.95);

Kfold = 5;
Kiter = N/Kfold;

D = zeros(lb*lK,1);
for k = 1:lK
    K1 = Kg1(k);
    K2 = Kg2(k);
    weight = get_weight_Naive(X,Si,K1,K2);
    weight=diag(weight);
    for j = 1:lb
        for i = 1:Kiter
            D((k-1)*lb+j,(i-1)*Kfold+1:i*Kfold) = diffi((i-1)*Kfold+1:i*Kfold,bw0(j),X,Y,Si,weight,q1t,q2t,q1x,q2x);
        end
    end
end 
G=mean(D,2);
idx = find(G==min(G), 1 );
idxj = mod(idx,lb);
if idxj == 0
    idxj=lb;
end
idxk = (idx-idxj)/lb +1;
bw = bw0(idxj);
K1 = Kg1(idxk);
K2 = Kg2(idxk);
end

function D = diffi(i,bw,X,Y,Si,weight,q1t,q2t,q1x,q2x)
%N = length(X);
Sic = Si(i,:);
Yc = Y(i,:);
Xc = X(i,:);
fc = zeros(1,length(Sic));
fc(Sic<=q2t & Sic>=q1t & min(Xc<=q2x,[],2) & min(Xc>=q1x,[],2))=1;


Si(i,:)=[];
Y(i,:)=[];

wc = weight(i);
weight(i)=[];

[Ethyes,~] = CaliEst_Naive(bw,Sic,Si,Y,weight);

D = (wc'.*Yc' - Ethyes).^2.*fc;

end

