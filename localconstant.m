function f = localconstant(x,X,Y)

N = length(X);
n = length(x);
x = reshape(x,1,n);
X = reshape(X,N,1);
Y = reshape(Y,N,1);

%Select bandwidth
w = X*0+1;
q1 = quantile(X,0.05);
q2 = quantile(X,0.95);
w(X<q1 | X>q2) = 0;

hROT = 1.06*std([X',Y'])*N^(-1/5);
hc = hROT/1.5:((hROT*2-hROT/1.5)/19):hROT*2;

ycv=arrayfun(@(i)lc(X,Y,hc(i)),1:length(hc),'UniformOutput',false);
ycv = cell2mat(ycv);
CV = sum((Y - ycv).^2.*w,1);
h = max(hc(CV==min(CV)));

%Fit the local constant estimator given bandwidth h
    d=x-X; %N x n
    z=normpdf(d./h);
    f=sum(z.*Y,1)./(sum(z,1));

end

function fcv = lc(X,Y,h)

N = length(X);
X = reshape(X,N,1);
Y = reshape(Y,N,1);

%Fit the local constant estimator given bandwidth h
    d=X'-X; %N x N
    z=normpdf(d./h);
   
    Dencv = sum(z,1) - z; %N x N
    Numcv = sum(z.*Y,1) - z.*Y; %N x N
    fcv = Dencv./Numcv;
    fcv = diag(fcv); %N x 1
    
end