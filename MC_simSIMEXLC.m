%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wei Huang and Zheng Zhang (2022).
% Nonparametric Estimation of the Continuous Treatment
% Effect with Measurement Error
%
% Main Matlab codes for Monte Carlo simulations on Estimation of ADRF
%
% *****************************************************
% Subroutines directly used in this main code:
%    (0) DGP1_ME.m; DGP2_ME.m; DGP11_ME.m; DGP12_ME.m: Data generation
%    processes
%    (1) Opt_Cali.m %theoretical optimal parameters
%    (2) Opt_Calibw.m %h0=hPI, K: by GCV and h: theoretical optimal
%    (3) SIMEXRegbw.m %our SIMEX method to choose h
%    (4) get_weightCon.m %estimate pi
%    (5) weightNWDecUknown.m %estimate mu
% *****************************************************
%
% Written by:
%    Wei Huang
%    Lecturer
%    School of Mathematics and Statistics, The University of Melbourne
%
% Last updated:
%    19 April, 2022.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear
%close all
%clc

tic;                                           % start timing
start_date = now;                              % start date
disp(['Job started: ', datestr(start_date)]);  % display start date

%% Step 1: Set-up
N=500; %sample size
J=200; %Number of Monte-Carlo samples
Bm=35; %Number of bootstrap in SIMEX method of selecting bandwidth

errortype = 'norm'; %Type of the error distribution, 2 choices: 'Lap' or 'norm'
        
DGP = 1;  %Data generation process: 1, 2, 11, 12

warning('off','all')

if DGP == 1
   [~,~, ~, ~,~, b,t0,mu,vUXr] = DGP1_ME(10000,1,errortype);
elseif DGP == 2
   [~,~, ~, ~,~, b,t0,mu,vUXr] = DGP2_ME(10000,1,errortype);
elseif DGP == 11
   [~,~, ~, ~,~, b,t0,mu,vUXr] = DGP11_ME(10000,1,errortype);
elseif DGP == 12
   [~,~, ~, ~,~, b,t0,mu,vUXr] = DGP12_ME(10000,1,errortype);
end

if strcmp(errortype,'Lap')==1
	varU = 2*b^2;
elseif strcmp(errortype,'norm')==1
	varU = b^2;
end
               
parfor j =1:J
    j
    %% Generate data
        if DGP == 1
            [Y,X, ~, Si,~, ~,~,~,~] = DGP1_ME(N,j,errortype);
        elseif DGP == 2
            [Y,X, ~, Si,~, ~,~,~,~] = DGP2_ME(N,j,errortype);
        elseif DGP == 11
            [Y,X, ~, Si,~, ~,~,~,~] = DGP11_ME(N,j,errortype);
        elseif DGP == 12
            [Y,X, ~, Si,~, ~,~,~,~] = DGP12_ME(N,j,errortype);
        end
    YDG(:,j) = Y;
    XDG(:,j) = X;
    SDG(:,j) = Si;
    
    %Y = YDG(:,j);
    %X = XDG(:,j);
    %Si = SDG(:,j);

    %Our method using theoretical optimal smoothing parameters (CM)
    [KOpt(j),bwOpt(j),hwOpt(j),ridgeOpt(j),idrOpt(j),idxOpt(j),ISEOpt(j)]...
        = Opt_Cali(t0,mu,Si,Y,X,errortype,b);
    weight = get_weightCon(t0,X,X,Si,errortype,hwOpt(j),KOpt(j),b);
    [EDOpt(j,:),~,~] = weightNWDecUknown(t0,Si,Y,weight,errortype,b,...
        bwOpt(j),ridgeOpt(j));
    
    %Choose h0 and K
    hPI(j) = PI_deconvUknownth4(Si,errortype,varU,b);
    Kg = 2:4;
    lhK = length(Kg);
    eX = exp(X);
    meX = mean(eX);
    GCV = zeros(3,1);
    for i = 1:lhK
        try
        weight = get_weightCon(t0,X,X,Si,errortype,hPI(j),Kg(i),b);
        %weight = zeros(N,length(t0));
        EXT = weightNWDecUknown(t0,Si,eX,weight,errortype,b,hPI(j),0.01);
        GCV(i) = sum((EXT - meX).^2)/(1-N^(-1)*Kg(i))^2;
        catch ME
        disp(ME.message)
        GCV(i) = Inf;
        end
    end
    K(j) = min(Kg(GCV==min(GCV)));
    weight = get_weightCon(t0,X,X,Si,errortype,hPI(j),K(j),b);
    
    %Our method using h0=hPI, K selected by method in section 5.2 and
    %theretical optimal h (CM_tilde)
    [bwOpttilde(j),~,ridgeOpttilde(j),~,~,ISEOptbw(j)]...
        = Opt_Calibw(t0,mu,Si,Y,X,errortype,b,K(j));
    [EDOptbw(j,:),~,~] = weightNWDecUknown(t0,Si,Y,weight,errortype,b,...
       bwOpttilde(j),ridgeOpttilde(j));
   
    
    %Our method using our data-driven smoothing parameters (CM_hat)
    [bwSIMEX(j),bwStar(j),bw1(:,j),bw2(:,j),ridgeSIMEX(j)] = SIMEXRegbw(Si,Y,X,K(j),errortype,b,Bm);
    [ENSIMEX(j,:),~,~] = weightNWDecUknown(t0,Si,Y,weight,errortype,b,bwSIMEX(j),ridgeSIMEX(j));
    ISESIMEX(j) = mean((ENSIMEX(j,:)-mu).^2);
end



filename=sprintf('DGP%d-%s-vUXr-%.2f-N%d-Opt-%s.mat',DGP,errortype,vUXr,N,date);
save(filename)

%% Step 3: Report computational time 
time = toc;        % finish timing
end_date = now;    % end date
disp('*****************************************');
disp(['Job started: ', datestr(start_date)]);
disp(['Job finished: ', datestr(end_date)]);
disp(['Computational time: ', num2str(time), ' seconds.']);
disp(['Computational time: ', num2str(time / 60), ' minutes.']);
disp(['Computational time: ', num2str(time / 3600), ' hours.']);
disp('*****************************************');
disp(' ');

load gong.mat;
sound(y)
