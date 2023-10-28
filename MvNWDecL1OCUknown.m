function CV=MvNWDecL1OCUknown(Ti,Si,X,Y,errortype,b,h,bw,rhogrid)

%Author: Wei Huang modified from Aurore Delaigle's code
	%Compute a version of weighted CV used in SIMEX, for binned data
	%Si plays the role of the contaminated data
	%midbin is the vector of the centers of the binned data that play the role of the non contaminated data

	%Default values of phiU(t)=characteristic function of the errors
	%If you want to consider another error type, simply replace phiU by the characteristic function of your error type

    dimX = dim(X,2);
    x = X';
    K = arrayfun(@(j)normpdf((x(j,:)-X(:,j))/bw(j)), 1:dimX, 'UniformOutput',false);
    
    Kx = 1;
    for j=1:dimX
        Kx = Kx.*cell2mat(K(j));
    end
    
	if strcmp(errortype,'Lap')==1
		phiU=@(t) 1./(1+b^2*t.^2);
	elseif strcmp(errortype,'norm')==1
		phiU = @(t) exp(-b^2*t.^2/2);
	end


	%phiK: Fourier transform of the kernel K. You can change this if you wish to use another kernel but make sure 
	%you change the range of t-values, which should correspond to the support of phiK		
	%The phiK used in this code must be the same as in the other codes (hSIMEX + NW codes)
	phiK = @(t) (1-t.^2).^3;
	dt = .001;
	t = (-1:dt:1)';
	th=t/h;
	longt=length(t);
	

	n=length(Si);

	
	%Compute the empirical characteristic function of W (times n) at t/h: \hat\phi_W(t/h)
	OO=outerop(Si,t/h,'*');
    csO=cos(OO);
    snO=sin(OO);
    clear OO;
	
	rehatphiW=Kx'*csO;
    imhatphiW=Kx'*snO;

	%Compute \sum_j Y_j e^{itW_j/h}
	Z = Y.*Kx;
    renum=Z'*csO;
    imnum=Z'*snO;

	%Compute \hat m(M_i) where M_i is the middle of the bin in which X_i (the non contaminated data) lies
	xt=outerop(t/h,Ti,'*');
    cxt=cos(xt);
    sxt=sin(xt);
    clear xt;
	

	phiUth=phiU(th);
	matphiKU=reshape(phiK(t)./phiUth,1,longt);
    
	Den=(rehatphiW.*matphiKU).*cxt'+(imhatphiW.*matphiKU).*sxt';
    Den = sum(Den,2);
    Den = Den';
	Num=(renum.*matphiKU).*cxt'+(imnum.*matphiKU).*sxt';
    Num = sum(Num,2);
    Num = Num';


	%Compute from there the leave-one-out version \hat m_{-i}(M_i)
	
	csO=csO';
	snO=snO';
	
	Den=Den-normpdf(0)*matphiKU*(csO.*cxt)-normpdf(0)*matphiKU*(snO.*sxt);
	for i=1:n
		csO(:,i)=csO(:,i)*Y(i);
		snO(:,i)=snO(:,i)*Y(i);
	end
	
	Num=Num-normpdf(0)*matphiKU*(csO.*cxt)-normpdf(0)*matphiKU*(snO.*sxt);


	%Finally compute weighted CV, where the ith term of the sum is weighted by f_W(W_i)
	bwp = prod(bw);
	rhogrid=rhogrid*(2*pi*h*bwp*n)/dt;
	%hW=1.06*sqrt(var(Si))*n^(-1/5);
	%xout=outerop(Si,Si,'-');
	%fWEF=normpdf(xout,0,hW)*ones(n,1)/n;
    %[fWEF,~,~]=ksdensity([X',W'],[X',W']);
    fWEF = zeros(n,1);
    q1 = quantile(Ti,0.05);
    q2 = quantile(Ti,0.95);
    fWEF(Ti<q2 & Ti>q1) = 1;
	
	CV=0*rhogrid;
	for krho=1:length(rhogrid)
		rho=rhogrid(krho);
		dd=Den;
		dd(dd<rho)=rho;

		mhatstar=Num./dd;
		mhatstar=reshape(mhatstar,1,n);		
		CV(krho)=sum(fWEF'.*(Y'-mhatstar).^2);
	end

