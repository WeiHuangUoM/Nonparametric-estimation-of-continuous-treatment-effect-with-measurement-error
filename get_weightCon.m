function weight = get_weightCon(t0,x,X,Si,errortype,h,K,b)
% Input: t0: n*1 vector of t values to evaluate E(Y*(t)) and pi_0(t,x)
%        x: nx*1 vector of x values to evaluate pi_0(t,x)
%        X: N*r matrix of the data of covariables X
%        Si: N*1 vector of the data of variable S
%        errortype: type of measrurement error -- Lap or norm
%        h: bandwidth h0 for estiamte pi_0(t,x)
%        K: number of seive basis functions to use
%        b: parameter to calculate f_U (density of measurement error)
% Output: weight: nx*n matrix of the estimators of pi(t,x)
% Sub-routines: kerU_d.m

N = length(Si);
ini = zeros(K,1);


vmat = arrayfun(@(j)FUNv(X(j)),1:N,'UniformOutput',false);
vmat = cell2mat(vmat');
vmat = vmat(:,1:K); %N*K matrix of v_{K}(X_i)'s

%orthogonomalize basis matrix
exp_vv = (1/N) * (vmat' * vmat);               % K2 x K2
exp_vv_half = chol(exp_vv, 'lower');
exp_vv_half_inv = exp_vv_half \ eye(K);
vmat_std = vmat * exp_vv_half_inv';

n = length(t0);
t0 = reshape(t0,1,n);
kerUvalue = kerU_d(t0,Si,errortype,b,h,0);
kerUvalue(kerUvalue<0) =0;
[LamHat,~] = arrayfun(@(i)optlam(ini,vmat_std,kerUvalue(:,i)),1:n,'UniformOutput',false);
LamHat = cell2mat(LamHat); %K x nt matrix of Lambda hat for t0

%vmat for output x
n = length(x);
vmat = arrayfun(@(j)FUNv(x(j)),1:n,'UniformOutput',false);
vmat = cell2mat(vmat');
vmat = vmat(:,1:K); %nx*K matrix of v_{K}(X_i)'s
vmat_std = vmat * exp_vv_half_inv';

term = vmat_std*LamHat;
weight = drho(term); %nx*nt matrix of rho'{LambdaHat_{t_j}^T*v_K(X_i)}
end


function [lamhat,fval] = optlam(ini,vmat,kerUvalue)
options= optimoptions('fminunc','Algorithm','quasi-newton','Display','off');
fun = @(lam)Obj_truncate(lam,vmat,kerUvalue);
[lamhat,fval] = fminunc(fun,ini,options);
end


function f = Obj_truncate(lam,vmat,kerUvalue)
normK = sum(kerUvalue);

meanvmat = mean(vmat,1);

K = length(lam);
Lam = reshape(lam,K,1);

u = Lam'*vmat';
Y = rho(u);

term1 = Y*kerUvalue;
term2 = Lam'*meanvmat';

f = -term1./normK + term2;

%if nargout > 1 % supply gradient
%    a = drho(u);
%    c = vmat'.*a;
%    term1 = c*kerUvalue;
%    term2 = meanvmat';
%    g = -term1./normK + term2;
%end

%if nargout > 2 % supply Hessian
%    a = d2rho(u);
%    cv = a.*kerUvalue'.*vmat';
%    term1 = cv*vmat;
%    H = -term1./normK;
%end

end

function r = rho(x)
r = -exp(-x-1);
end

function dr = drho(x)
dr = exp(-x-1);
end


function bv = FUNv(x)

bv = [x.^0,x,x.^2,x.^3,x.^4,x.^5,x.^6,x.^7,x.^8,x.^9,x.^10];
end

