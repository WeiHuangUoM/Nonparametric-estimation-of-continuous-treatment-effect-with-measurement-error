function weight = get_weightCon(t0,x,X,Si,errortype,h,KX,b)
% Input: t0: n*1 vector of t values to evaluate E(Y*(t))
%        ini: initial value for maximising G to find Lambda Hat
%        X: N*r matrix of the data of covariables X
%        Si: N*1 vector of the data of variable S
%        errortype: type of measrurement error -- Lap or norm
%        h: bandwidth for calculating K_U
%        K: number of seive basis functions to use
%        b: parameter to calculate f_U (density of measurement error)
% Output: weight: N*n matrix of the estimators of pi(t_j,X_i)

N = length(Si);
dimX = size(X,2);
K = 1 + (KX-1)*dimX;
ini = zeros(1,K);

vmat = arrayfun(@(j)FUNv(X(j,:)),1:N,'UniformOutput',false);
vmat = cell2mat(vmat');
vmat = vmat(:,1:K); %N*K matrix of v_{K}(X_i)'s

n = length(t0);
t0 = reshape(t0,1,n);
kerUvalue = kerU(t0,Si,errortype,b,h);

[LamHat,~] = arrayfun(@(i)optlam(ini,vmat,kerUvalue(:,i)),1:n,'UniformOutput',false);
LamHat = cell2mat(LamHat'); %nt*K matrix of Lambda hat for t0

%vmat for output x
n = size(x,1);
vmat = arrayfun(@(j)FUNv(x(j,:)),1:n,'UniformOutput',false);
vmat = cell2mat(vmat');
vmat = vmat(:,1:K); %nx*K matrix of v_{K}(X_i)'s

term = vmat*LamHat';
weight = drho(term); %nx*nt matrix of rho'{LambdaHat_{t_j}^T*v_K(X_i)}
end


function [lamhat,fval] = optlam(ini,vmat,kerUvalue)

options= optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'HessianFcn','objective','Display','off');
%options= optimoptions('fminunc','Algorithm','quasi-newton','Display','off');
fun = @(lam)Obj(lam,vmat,kerUvalue);
[lamhat,fval] = fminunc(fun,ini,options);


end

function [f,g,H] = Obj(lam,vmat,kerUvalue)

normK = sum(kerUvalue);
meanvmat = mean(vmat,1);

K = length(lam);
Lam = reshape(lam,K,1);

u = Lam'*vmat';
Y = rho(u);

term1 = Y*kerUvalue;
term2 = Lam'*meanvmat';

f = -term1./normK + term2;

if nargout > 1 % supply gradient
    a = drho(u);
    c = vmat'.*a;
    term1 = c*kerUvalue;
    term2 = meanvmat';
    g = -term1./normK + term2;
end

if nargout > 2 % supply Hessian
    a = d2rho(u);
    cv = a.*kerUvalue'.*vmat';
    term1 = cv*vmat;
    H = -term1./normK;
end

end


function r = rho(x)
r = -exp(-x-1);
end

function dr = drho(x)
dr = exp(-x-1);
end

function d2r = d2rho(x)
d2r = -exp(-x-1);
end

function bv = FUNv(x)

bv = [1,x];
end

