function [f,KerZ,LU] = fdecCondUKnown(N,t0,Si,Z,h1,h2,errortype,sigU,l)

%Compute (the derivative of) the deconvolution kernel of T=t0 given a
%multivariate covariate Z=Z

%N: sample size
%h1: bandwidth for Si
%h2: bandwidth for Z
%l: order of derivative, e.g. l=1: first derivative
%
%errortype: 'Lap' for Laplace errors and 'norm' for normal errors. For other error distributions, simply redefine phiU below 
%sigU: parameter of Laplace or normal errors used only to define phiU.
%
%Sub-routines: (1) kerU_d.m

if isempty(Z) 
    KerZ = ((1:N).^0)'; %N x 1 vector of 1.
else
%Reshape
Z = reshape(Z,N,[]);
dz = size(Z,2);
z = reshape(Z,[],dz);
nz = size(z,1);
t0 = reshape(t0,1,length(t0));
Si = reshape(Si,N,1);

%measure the L2 distance between Z and z
d = zeros(N,nz);
for j = 1:nz
    d(:,j) = sqrt(sum((Z - z(j,:)).^2,2));
end
%normal kernel of Z
KerZ=normpdf(d./h2)/h2; %N x nz
end

LU = kerU_d(t0,Si,errortype,sigU,h1,l)*N; %N x lt

Num = KerZ'*LU; %nz x lt
Den = (sum(KerZ,1))'; %nz x 1

f = Num./Den; %nz x lt

if l == 10
%Normalization
LU0 =  kerU_d(t0,Si,errortype,sigU,h1,0)*N;
Num0 = KerZ'*LU0; %nz x lt

f0 = Num0./Den; %nz x lt
scale = sum(f0,2)*(t0(2)-t0(1));

f = f./scale;
end

end