function weight = get_weight_Naive(X,Si,K1,KX)
% Input: ini: initial value for maximising G to find Lambda Hat
%        X: N*r vector of the data of variable X
%        Si: N*1 vector of the data of variable S
% Output: weight: N*1 vector of the estimators of pi_i's

DimX = size(X,2);
K2 = 1+(KX-1)*DimX;

ini = zeros(K1,K2);

N=length(Si);

umat = arrayfun(@(j)FUNu(Si(j)),1:N,'UniformOutput',false);
umat = cell2mat(umat'); %N*K_1 matrix of the u_{K_1}(S_i)'s
umat = umat(:,1:K1);
vmat = arrayfun(@(j)FUNv(X(j,:)),1:N,'UniformOutput',false);
vmat = cell2mat(vmat'); %N*K_2 matrix of the v_{K_2}(X_i)'s
vmat = vmat(:,1:K2);

[LamHat,~] = optlam(ini,umat,vmat,N); %Output the (K_1*K_2)*1 vector of the
                                      %estimator of Lambda 
%LamHat = reshape(LamHat,K1,K2);       %reshape Lambda Hat to a K_1*K_2 matrix.

u = umat*LamHat*vmat'; 
v = u';  %N*1 vector of u_{K_1}(S_i)^T*LambdaHat*v_{K_2}(X_i)

weight = drho(v); %obtain the weight pi_i
end



function [lamhat,fval] = optlam(ini,umat,vmat,N)

options= optimoptions('fminunc','Algorithm','trust-region','StepTolerance',1e-10,'SpecifyObjectiveGradient',true,'Display','off');
%options= optimoptions('fminunc','Algorithm','quasi-newton','Display','off');
fun = @(lam)Obj(lam,umat,vmat,N);
[lamhat,fval] = fminunc(fun,ini,options);


end


function [f,g] = Obj(Lam,umat,vmat,N)

meanumat = mean(umat,1);
meanvmat = mean(vmat,1);

%Lam = reshape(lam,K1,K2);

u = umat*Lam.*vmat;
v = sum(u,2);
Y = rho(v);

term1 = mean(Y);
term2 =meanumat*Lam*meanvmat';

f = -term1 + term2; %Function G(Lambda)

if nargout > 1 % supply gradient
    a = drho(v);
    c = vmat.*a;
    term1 = umat'*c/N;
    term2 = meanumat'*meanvmat;
    g = -term1 + term2; %Gradient of G(Lambda)
end

end

function r = rho(x)
r = -exp(-x-1);
end

function dr = drho(x)
dr = exp(-x-1);
end


function bv = FUNu(x)

bv = [1,x];
end

function bv = FUNv(x)

bv = [1,x];
end




