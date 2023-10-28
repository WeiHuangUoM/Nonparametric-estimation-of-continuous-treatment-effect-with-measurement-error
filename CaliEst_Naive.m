function [Ethyes,den]= CaliEst_Naive(bw,t0,Si,Y,weight)

t0 = reshape(t0,1,length(t0));

kt = (t0-Si)/bw; %N*n
kervalue = normpdf(kt); 
sumK = sum(kervalue,1);

num = sum(weight.*Y.*kervalue,1)./sumK; %1*n
den = sum(weight.*kervalue,1)./sumK;

Ethyes = num./(den); %1*n

end
