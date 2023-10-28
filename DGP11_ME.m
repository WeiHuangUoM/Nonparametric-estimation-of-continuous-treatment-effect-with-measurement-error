function [Y,X, Ti, Si,Uy, b,t0,mu,vUXr] = DGP11_ME(N,Seed,errortype)

rng(Seed);
Uy = normrnd(0,1,[N,1]);
sdT = sqrt(1.016);

X = 0.2+rand(N,1)*0.6;
Ut = normrnd(0,1,[N,1]);

Ti = sqrt(X) -0.7 + Ut;

if strcmp(errortype,'Lap')==1
    vUXr = 0.25;
	b = sqrt(vUXr*sdT^2)/sqrt(2);
    Si = Ti + (rlap(b,1,N))';
elseif strcmp(errortype,'norm')==1
    vUXr = 0.2;
	b = sqrt(vUXr*sdT^2);
    Si = Ti + normrnd(0,b,[N,1]);
end

Y = Ti + exp(X) +Uy;

mint = quantile(Ti,0.1);
maxt = quantile(Ti,0.9);
t0 = mint:(maxt-mint)/99:maxt;
mu = t0 + mean(exp(X) +Uy);

end