function [Y,X, Ti, Si,Uy, b,t0,mu,vUXr] = DGP2_ME(N,Seed,errortype)

rng(Seed);
Uy = rand(N,1);
sdT = sqrt(0.18/12+1);

X = sum(0.3*rand(N,2),2);
Ut = normrnd(0,1,[N,1]);

Ti = X.^2 + 1  + Ut;

if strcmp(errortype,'Lap')==1
    vUXr = 0.25;
	b = sqrt(vUXr*sdT^2)/sqrt(2);
    Si = Ti + (rlap(b,1,N))';
elseif strcmp(errortype,'norm')==1
    vUXr = 0.2;
	b = sqrt(vUXr*sdT^2);
    Si = Ti + normrnd(0,b,[N,1]);
end

Y = exp(-6+6*Ti)./(1+exp(-6+6*Ti)) + X +Uy;

mint = quantile(Ti,0.1);
maxt = quantile(Ti,0.9);
t0 = mint:(maxt-mint)/99:maxt;
mu =  exp(-6+6*t0)./(1+exp(-6+6*t0)) + mean(X +Uy);

end