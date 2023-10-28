%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wei Huang and Zheng Zhang (2021).
%  Wei Huang and Zheng Zhang (2022).
% Nonparametric Estimation of the Continuous Treatment
% Effect with Measurement Error
%
% Main Matlab codes for Monte Carlo simulations on Estimation of ADRF
% using Naive estimator.
% *****************************************************
% Subroutines directly used in this main code:
%    (0) DGP1_ME.m; DGP2_ME.m; DGP11_ME.m; DGP12_ME.m: Data generation
%    processes
%    (1) Opt_Naive.m %theoretical optimal parameters
%    (2) CV_NaiveFull.m %10-fold cross-validation for parameters
%    (3) get_weight_Naive.m %naive way to estimate pi
%    (4) CaliEst_Naive.m %naive way to estimate mu
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
% close all
%clc

tic;                                           % start timing
start_date = now;                              % start date
disp(['Job started: ', datestr(start_date)]);  % display start date

%% Step 1: Set-up
N=750;
J=200;

errortype = 'Lap'; %Type of the error distribution, 2 choices: 'Lap' or 'norm'
        
DGP = 1;

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

    %Naive estimator using theoretical optimal smoothing parameters (NV)
    [bwNOpt(j),idxNOpt(j),K1NOpt(j),K2NOpt(j),ISENOpt(j)] = Opt_Naive(t0,mu,Si,Y,X);
    weight = get_weight_Naive(Si,X,X,Si,K1NOpt(j),K2NOpt(j));
    weight=diag(weight);
    [NV(j,:),~]= CaliEst_Naive(bwNOpt(j),t0,Si,Y,weight);
    
    %CV K1 K2 bw for Naive
    [bwFNCV(j),K1NCV(j),K2NCV(j),~] = CV_NaiveFull(Si,Y,X);
    weight = get_weight_Naive(Si,X,X,Si,K1NCV(j),K2NCV(j));
    weight=diag(weight);
    ENFCV(j,:)= CaliEst_Naive(bwFNCV(j),t0,Si,Y,weight);
    ISENFCV(j) = mean((ENFCV(j,:)-mu).^2);
    
  
end

filename=sprintf('DGP%d-%s-vUXr-%.2f-N%d-Naive-%s.mat',DGP,errortype,vUXr,N,date);
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

%load gong.mat;
%sound(y)

