% This demo shows how to use DivSBL to reconstruct real-world audio signal.

close all
clear; clc;
%Set your own path
addpath(genpath('C:\Users\king\Desktop\DivSBL_public'));


% info =audioinfo('test_24k.wav');
%[audio,Fs] = audioread('test_24k.wav');
info =audioinfo('-0SdAVK79lg.wav');
[audio,Fs] = audioread('-0SdAVK79lg.wav');

audiolength = 480;
t = 1:1:audiolength;
audio = audio(:,1);
audio = audio(t);
figure(1),
plot(t,audio);
xlabel('Time');
ylabel('Audio Signal');

rng(0,'twister') 
M = 150; 
N = audiolength;
SNR = 15;  

Wgen = dct(audio);
iterNum = 10;    % number of experiments


for it = 1: iterNum

h_em_l = 16;
fprintf('\n\nRunning %d:\n',it);

% Generate the known matrix with columns draw uniformly from the surface of a unit hypersphere
Phi = randn(M,N);
Phi = Phi./(ones(M,1)*sqrt(sum(Phi.^2)));


%% Measurements
% noiseless signal
signal = Phi * Wgen;

% Observation noise   
stdnoise = std(signal)*10^(-SNR/20);
noise = randn(M,1) * stdnoise;

% Noisy signal
Y = signal + noise;

IND = find(abs(Wgen)>0);


% ======================================================================
%             Algorithm Comparison
% ======================================================================


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Benchmark
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
supt = find( abs(Wgen)>0);
x_ls = pinv(Phi(:,supt)) * Y;
mse_bench(it) = (norm(Wgen(supt) - x_ls,'fro')/norm(Wgen,'fro'))^2; 
fprintf('   bench:  MSE: %g;\n', mean(mse_bench));



%% ================== BSBL-EM ==============
blkStartLoc = [1:h_em_l:N];
learnlambda = 1;         % when SNR < 25dB, set learnlambda = 1;
                         % when SNR >=25dB, set learnlambda = 2;
                         % when no noise present, set learnlambda = 0;
tic;

Result1 = BSBL_EM(Phi,Y,blkStartLoc,learnlambda);

t_em(it) = toc;
mse_em(it) = (norm(Wgen - Result1.x(:,Result1.count),'fro')/norm(Wgen,'fro'))^2;


fprintf('BSBL-EM(learn correlation) : time: %4.3f, MSE: %g\n',mean(t_em),mean(mse_em));


%% ================== DivSBL ==============
blkStartLoc_l = [1:h_em_l:N];
tic;
Result_DV = DivSBL(Phi, Y, blkStartLoc_l);

t_em_DV(it) = toc;
mse_em_DV(it) = (norm(Wgen - Result_DV.x(:,Result_DV.count),'fro')/norm(Wgen,'fro'))^2;


fprintf('DivSBL: time: %4.3f, MSE: %g\n',mean(t_em_DV),mean(mse_em_DV));


%% ================== PC-SBL ==============

% Initialization of parameters  
n=N;                                          % signal dimension
m=M;                                           % number of measurements
A=Phi;
y=Y;

a=0.5;
b=1e-10;
c=1e-10;
d=1e-10;
eta=1;

tic;
iter=0;
iter_mx=500;
D=eye(n);
sigma2=1;
alpha_new=ones(n,1);
var_new=inv(A'*A/sigma2+D);
mu_old=ones(n,1);
mu_new=1/sigma2*var_new*A'*y;
gamma_new=1/sigma2;
while iter<iter_mx& norm(mu_new-mu_old)>1e-6
iter=iter+1;
mu_old=mu_new;
mul=[mu_new(2:n);0];
mur=[0;mu_new(1:n-1)];
var=diag(var_new);
varl=[var(2:n);0];
varr=[0;var(1:n-1)];
E=mu_new.^2+eta*mul.^2+eta*mur.^2+var+eta*varl+eta*varr;
alpha_new=a./(0.5*E+b);
idx1=find(alpha_new>1e10);
alpha_new(idx1)=1e10;
alf=[alpha_new(2:n); 0];                                %   left-shifted version
arf=[0; alpha_new(1:n-1)];                              %   right-shifted version
D=diag(alpha_new+eta*alf+eta*arf);
%=============================================
%  estimate the variance of noise
 num=(y-A*mu_old)'*(y-A*mu_old)+trace(var_new*A'*A)+2*d;
den=m+2*c;
sigma2=num/den;
%==============================================
var_new=inv(A'*A/sigma2+D);
mu_new=1/sigma2*var_new*A'*y;
end
x_new=mu_new;
t_PC(it)=toc;
mse_PC(it)=(norm(Wgen - x_new,'fro')/norm(Wgen,'fro'))^2;
fprintf('PC-SBL: time: %4.3f, MSE: %g\n',mean(t_PC),mean(mse_PC));


%% ==================SBL(RVM) ==============
% solve by BCS
initsigma2 = std(Y)^2/1e2;
tic;
[weights,used,sigma2,errbars] = BCS_fast_rvm(A,Y,initsigma2,1e-8);
t_BCS(it) = toc;
%fprintf(1,'BCS number of nonzero weights: %d\n',length(used));
x_BCS = zeros(N,1); err_BCS = zeros(N,1);
x_BCS(used) = weights; err_BCS(used) = errbars;

mse_SBL(it) = (norm(Wgen - x_BCS,'fro')/norm(Wgen,'fro'))^2;


fprintf('SBL: time: %4.3f, MSE: %g\n',mean(t_BCS),mean(mse_SBL));

%% ================== Group Lasso ==============
tic;
cvx_begin quiet
variable x_group_lasso(N,1)
k_cvx=h_em_l; 
mu=3*1e-1;
sum_of_norms = 0;

for i = 1:k_cvx:(N)
indices = i:min(i+k_cvx-1, N); 
x_subset = x_group_lasso(indices);
norm_value = norm(x_subset, 2);
sum_of_norms = sum_of_norms + norm_value;
end
minimize(square_pos(norm(A * x_group_lasso - Y, 'fro')) + mu* sum_of_norms);
cvx_end
t_group_lasso(it)=toc;

mse_group_lasso(it) = (norm(Wgen - x_group_lasso,'fro')/norm(Wgen,'fro'))^2;


fprintf('Group Lasso(CVX): time: %4.3f, MSE: %g\n',mean(t_group_lasso),mean(mse_group_lasso));

%% ================== Group BPDN ==============

sigma = 1e-6;
tic;
numElements = N / h_em_l;
repeatedArray = repmat(1:numElements, h_em_l, 1);
groups = repeatedArray(:);
x_GBPDN = spg_group( A, Y, groups, sigma);
t_GBPDN(it) = toc;

mse_GBPDN(it) = (norm(Wgen - x_GBPDN,'fro')/norm(Wgen,'fro'))^2;


fprintf('Group BPDN: time: %4.3f, MSE: %g\n',mean(t_GBPDN),mean(mse_GBPDN));





%% ================== StructOMP ==============

H=N; W=1; n=H*W; P=1; 
K=100; %total non-zero num %Pre-set
q=3;% block Num  
Edge4=Edge4Index(H, W);
if W==1
    Edge4=Edge4(:,1:3);
end

[B, Bm]=GetBlocksMatrix(H, W, 1);
[BC, BCm]=GetBlocksConnectionMatrix(Bm, Edge4);

mq=length(B);
 tic;
    lamada=1;
    cl0=2*lamada*q*log2(mq)+K;
    [x_StrOMP, input, norm_save6] = GraphOMP_CS(cl0,Phi,y, Bm, BCm, lamada, []);
    t_StrOMP=toc;

mse_StrOMP(it) = (norm(Wgen - x_StrOMP,'fro')/norm(Wgen,'fro'))^2;
fprintf('StrOMP: time: %4.3f, MSE: %g\n',mean(t_StrOMP),mean(mse_StrOMP));


end






%% Figures
figure(1);%Original

stem(Wgen, 'LineWidth', 1.7, 'Marker', 'none', 'Color', '#0072BD'); 
xlim([min(Wgen), N]);

title('Original', 'FontSize', 14);


%%
figure(2);%DivSBL

stem(Result_DV.x(:,Result_DV.count), 'LineWidth', 1.7, 'Marker', 'none', 'Color', '#0072BD'); 
xlim([min(Result_DV.x(:,Result_DV.count)), N]);

title('DivSBL', 'FontSize', 14); 

%%
figure(3);%BSBL-EM

stem(Result1.x(:,Result1.count), 'LineWidth', 1.7, 'Marker', 'none', 'Color', '#0072BD');
xlim([min(Result1.x(:,Result1.count)), N]);

title('BSBL-EM', 'FontSize', 14); 

%% 
figure(4);%PC-SBL

stem(x_new, 'LineWidth', 1.7, 'Marker', 'none', 'Color', '#0072BD'); 
xlim([min(x_new), N]);

title('PC-SBL', 'FontSize', 14); 

%%
figure(5);%SBL

stem(x_BCS, 'LineWidth', 1.7, 'Marker', 'none', 'Color', '#0072BD'); 
xlim([min(x_BCS), N]);

title('SBL', 'FontSize', 14); 


%%
figure(6);%Group Lasso

stem(x_group_lasso, 'LineWidth', 1.7, 'Marker', 'none', 'Color', '#0072BD'); 
xlim([min(x_group_lasso), N]);

title('Group Lasso', 'FontSize', 14); 


%%
figure(7);% Group BPDN

stem(x_GBPDN, 'LineWidth', 1.7, 'Marker', 'none', 'Color', '#0072BD');
xlim([min(x_GBPDN), N]);

title('Group BPDN', 'FontSize', 14); 

%%
figure(8);% StructOMP

stem(x_StrOMP, 'LineWidth', 1.7, 'Marker', 'none', 'Color', '#0072BD'); 
xlim([min(x_BOMP), N]);

title('StructOMP', 'FontSize', 14); 


