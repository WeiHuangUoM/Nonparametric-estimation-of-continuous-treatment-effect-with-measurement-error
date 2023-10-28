function reponse=rlap(szC,n1,n2)
%Author: Aurore Delaigle
% Generate a matrix of size n1 x n2 from a Laplace(szC)
y=unifrnd(0,1,n1,n2);
reponse= szC*log(2*y);
reponse(y>0.5)=-szC*log(2-2*y(y>0.5));

end

