% Recovery of block sparse images using DivSBL.

clear;close all;
%Set your own path here
addpath(genpath('C:\Users\king\Desktop\DivSBL_public'));
rng(5,'twister') 
N=128;
%parrots image (Set your own path here)
load('C:\Users\king\Desktop\DivSBL_public\Image demo\CS_test_images\Parrots.mat');
img=cdata;

ss=[N,N]; % size of image
X0=imresize(img,ss);
X0=double(X0);
X=X0;
[a,b]=size(X);

% Discrete Wavelet Transform
load DWTM.mat  
M=64;
R=randn(M,a);
SNR=25;
measure=R*X;
% Observation noise
stdnoise = std(measure)*10^(-SNR/20);
noise = randn(M,1) * stdnoise;

Y=measure+noise;

Phi=R*ww';

% figure(1)
% X=reshape(X,ss);
% imshow(uint8(X));
% title('Original Image')
% 
% figure(2)
% sp1 = inv(ww')*X;
% imshow(uint8(sp1));
% title('Sparse Image')
% 
% figure(3)
% plot(sp1(:,[1:5]))

count = 1;
%% =========================================================
 %                Proposed PC-SBL algorithm
 %=========================================================
X3=zeros(a,b);
eta=1;
tic;
for i=1:N
    rec1=PCSBL(Y(:,i),Phi,2,eta); % recover the image column by column
    X2(:,i)=rec1;
end
t_PC = toc;
X2=ww'*X2;  %  inverse-DWT transform
X22=reshape(X2,ss);
ERR=sqrt(sum(sum((X22-X0).^2,1),2)/sum(sum(X.^2,1),2));
fprintf('PC-SBL: time: %4.3f, NMSE: %g\n',mean(t_PC),mean(ERR));
figure(5)
imshow(uint8(X22));
title('PC-SBL');



%% ================== BSBL-EM  ==============

    h_em = 4;
    blkStartLoc = [1:h_em:N];
    learnlambda = 1;         % when SNR < 25dB, set learnlambda = 1;
                             % when SNR >=25dB, set learnlambda = 2;
                             % when no noise present, set learnlambda = 0;
    tic;
    
    for i=1:N
        Result1=BSBL_EM(Phi,Y(:,i),blkStartLoc,learnlambda); 
        X_BSBL(:,i)=Result1.x(:,Result1.count);
    end
    
    t_em= toc;
    X_BSBL=ww'*X_BSBL;  %  inverse-DWT transform
    X_BSBL=reshape(X_BSBL,ss);
    mse_em = sqrt(sum(sum((X_BSBL-X0).^2,1),2)/sum(sum(X.^2,1),2));
    
    
    fprintf('BSBL-EM(learn correlation) : time: %4.3f, NMSE: %g\n',mean(t_em),mean(mse_em));
    figure(6)
    imshow(uint8(X_BSBL));
    title('BSBL');


%%       =========================================================
%                    DivSBL
%     =========================================================
    h_em_l = 8;
    blkStartLoc_l = [1:h_em_l:N];
    tic;
    for i=1:N
    Result_DV = DivSBL(Phi,Y(:,i),blkStartLoc_l); 
    X_Div(:,i) = Result_DV.x(:,Result_DV.count);
    end
    
    t_em_DV= toc;
    X_Div=ww'*X_Div;  %  inverse-DWT transform
    X_Div=reshape(X_Div,ss);
    mse_em_DV = sqrt(sum(sum((X_Div-X0).^2,1),2)/sum(sum(X.^2,1),2));
   

    fprintf('DivSBL: time: %4.3f, MSE: %g\n',mean(t_em_DV),mean(mse_em_DV));

    figure(7)
    imshow(uint8(X_Div));
    title('DivSBL');

%% SBL(RVM)
% solve by BCS
initsigma2 = std(Y).^2./1e2;
tic;
 for i=1:N
    [weights,used,sigma2,errbars] = BCS_fast_rvm(Phi,Y(:,i),initsigma2(i),1e-8); 
    X_BCS = zeros(N,1); err_BCS = zeros(N,1);
    X_BCS(used) = weights; err_BCS(used) = errbars;
    X_SBL(:,i) = X_BCS;
 end

    t_BCS = toc;
    X_SBL=ww'*X_SBL;  %  inverse-DWT transform
    X_SBL=reshape(X_SBL,ss);

mse_SBL= sqrt(sum(sum((X_SBL-X0).^2,1),2)/sum(sum(X.^2,1),2));

fprintf('SBL: time: %4.3f, MSE: %g\n',mean(t_BCS),mean(mse_SBL));

    figure(8)
    imshow(uint8(X_SBL));
    title('SBL');

%% Group Lasso
tic;
for j=1:N
cvx_begin quiet
variable x_group_lasso(N,1)
k_cvx=h_em_l; 
    mu=3*1e-1;
sum_of_norms = 0;

for i = 1:k_cvx:(a)
    indices = i:min(i+k_cvx-1, a);
    x_subset = x_group_lasso(indices);
    norm_value = norm(x_subset, 2);
    sum_of_norms = sum_of_norms + norm_value;
end
minimize(square_pos(norm(Phi * x_group_lasso - Y(:,j), 'fro')) + mu* sum_of_norms);
cvx_end
X_group_lasso(:,j) = x_group_lasso;
end
t_group_lasso=toc;
    X_group_lasso=ww'*X_group_lasso;  %  inverse-DWT transform
    X_group_lasso=reshape(X_group_lasso,ss);

mse_group_lasso= sqrt(sum(sum((X_group_lasso-X0).^2,1),2)/sum(sum(X.^2,1),2));


fprintf('Group Lasso(CVX): time: %4.3f, MSE: %g\n',mean(t_group_lasso),mean(mse_group_lasso));

    figure(9)
    imshow(uint8(X_group_lasso));
    title('Group Lasso');

%% Group BPDN
sigma = 1e-6;
for i=1:N
tic;
numElements = a / h_em_l;
repeatedArray = repmat(1:numElements, h_em_l, 1);
groups = repeatedArray(:);
x_GBPDN = spg_group( Phi, Y(:,i), groups, sigma);
X_GBPDN(:,i) = x_GBPDN;
end
t_GBPDN = toc;
   X_GBPDN=ww'*X_GBPDN;  %  inverse-DWT transform
   X_GBPDN=reshape(X_GBPDN,ss);

mse_GBPDN= sqrt(sum(sum((X_GBPDN-X0).^2,1),2)/sum(sum(X.^2,1),2));


fprintf('Group BPDN: time: %4.3f, MSE: %g\n',mean(t_GBPDN),mean(mse_GBPDN));

    figure(10)
    imshow(uint8(X_GBPDN));
    title('Group BPDN');

%% StructOMP
H=a; W=1; n=H*W; P=1; 
K=50; %total non-zero num %Pre-set
q=3;% block Num  
Edge4=Edge4Index(H, W);
if W==1
    Edge4=Edge4(:,1:3);
end

[B, Bm]=GetBlocksMatrix(H, W, 1);
[BC, BCm]=GetBlocksConnectionMatrix(Bm, Edge4);

mq=length(B);
 tic;
 for i=1:N
    lamada=1;
    cl0=2*lamada*q*log2(mq)+K;
    [x_StrOMP, input, norm_save6] = GraphOMP_CS(cl0,Phi,Y(:,i), Bm, BCm, lamada, []);
    X_StrOMP(:,i) = x_StrOMP;
 end
    t_StrOMP=toc;
   X_StrOMP=ww'*X_StrOMP;  %  inverse-DWT transform
   X_StrOMP=reshape(X_StrOMP,ss);

mse_StrOMP= sqrt(sum(sum((X_StrOMP-X0).^2,1),2)/sum(sum(X.^2,1),2));

fprintf('StructOMP: time: %4.3f, MSE: %g\n',mean(t_StrOMP),mean(mse_StrOMP));

    figure(10)
    imshow(uint8(X_StrOMP));
    title('StructOMP');
