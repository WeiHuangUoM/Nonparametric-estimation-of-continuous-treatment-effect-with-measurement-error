meP = 3/4; %percentage in observation due to measurement error 
            %(1/6; 3/7; 3/4)

%Data preparation
load('HanebookData')
Y = Data(:,1);
Age = Data(:,2);
PIR = log(Data(:,3));
BMI = Data(:,4);
ALC = Data(:,5);
FH = Data(:,6);
AM = Data(:,7)-1;
PM = Data(:,8)-1;
Race = Data(:,9)-1;
Cal = Data(:,11);
X = [Age,PIR,BMI,Cal,ALC,FH,AM,PM,Race];
Si = log(5+Data(:,10)*100);

errortype='norm'; %Set errortype
varS = var(Si);

varU = meP*varS; %variance of measurement error
if strcmp(errortype,'Lap')==1
	b = sqrt(varU/2);
elseif strcmp(errortype,'norm')==1
	b = sqrt(varU);
end

Bm=35; %Set number of simulated data set in SIMEX
K=2; %Set K for calibration deconvolution method

%Set plot grid
q1 = quantile(Si,0.1);
q2 = quantile(Si,0.9);
t0 = q1:(q2-q1)/49:q2;



%Naive Estimator
[bwN,K1N,K2N,~] = CV_Naive(Si,Y,X);
weightN = get_weight_Naive(X,Si,K1N,K2N);
weightN=diag(weightN);
EN= CaliEst_Naive(bwN,t0,Si,Y,weightN);

%Deconvolution-Calibrationn Estimator
hPI = PI_deconvUknownth4(Si,errortype,varU,b);
weight = get_weightCon(t0,X,X,Si,errortype,hPI,K,b);
[bwD,bw1,bw2,ridgeD,bwStar] = SIMEXRegbw(Si,Y,X,errortype,b,Bm,K);
[ED,~] = weightNWDecUknown(t0,Si,Y,weight,errortype,b,bwD,ridgeD);


%filename = sprintf('NHANES_D_meP%.2f_T95_Cal.mat',meP);
%save(filename)
