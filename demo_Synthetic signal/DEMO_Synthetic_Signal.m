% This demo shows how to use DivSBL to reconstruct the synthetic signal.

close all
clear; clc;
%Set your own path
addpath(genpath('C:\Users\king\Desktop\DivSBL_public'));
rng(0,'twister')
M = 80;          % row number of the dictionary matrix 
N = 162;          % column number
blkNum = 4;       % nonzero block number
blkLen = 10;       % block length

SNR = 15;         % Signal-to-noise ratio
iterNum = 20;    % number of experiments

for it = 1: iterNum
    fprintf('\n\nRunning %d:\n',it);
    h_em_l = 18;       % harmless pre-set block length: N must divisible by it. 
    % Check if N is divisible by h_em_l
    if mod(N, h_em_l) ~= 0
        fprintf('Please enter an h_em_l that can divide N\n');
    end
    % Generate the known matrix with columns draw uniformly from the surface of a unit hypersphere
    Phi = randn(M,N);
    Phi = Phi./(ones(M,1)*sqrt(sum(Phi.^2)));


%% generate nonzero block coefficients

%% %    'Homoscedastic' Data Given by [1]
%     beta = ones(blkNum,1)*0.99;
%     blks(:,1) = randn(blkNum,1);
%     for i = 2 : blkLen
%         blks(:,i) = beta .* blks(:,i-1) + sqrt(1-beta.^2).*(ones(blkNum,1).*randn(blkNum,1));
%     end
%      % normalize along row and vectorize
%     nonzeroCoeff = reshape(blks',blkNum*blkLen,1);

%     %========================================================================
%     % put blocks at random locations and align with block partition (no overlap) 
%     %========================================================================
%     psbBlockNum = floor(N/blkLen);
%     ind2 = randperm(psbBlockNum);
%     indice2 = ind2(1:blkNum);
%     Wgen = zeros(N,1);  blkLoc2 = [];
%     for i = 1 : blkNum
%         Wgen( (indice2(i)-1)*blkLen + 1 : indice2(i)*blkLen ) = nonzeroCoeff((i-1)*blkLen + 1 : i*blkLen);
%         blkLoc2 = [blkLoc2,(indice2(i)-1)*blkLen + 1 : indice2(i)*blkLen];
%     end


%%   Generate 'Heteroscedastic' data
    
    % Generate each block size randomly.
    minGap = 2; 
    partitionIndices = sort(randperm(blkLen*blkNum-1-(blkNum-1)*minGap, blkNum-1));
    partitionIndices = partitionIndices + (1:blkNum-1)*minGap;
    
    resultArray = diff([0 partitionIndices blkLen*blkNum]);
  
    if resultArray(end) <= minGap-1
        resultArray(end-1) = resultArray(end-1) - minGap+1;
        resultArray(end) = resultArray(end)+1;  
    end
    
    nums = resultArray;


    % Randomly select blkNum positions to place the blocks generate above.
    ind2 = zeros(1, blkNum);
    for i = 1:blkNum
        ind2(i) = randi([1, N-max(nums)+1]);
    end
    ind2 = sort(ind2);
    for i = 1:blkNum-1 
        sub = ind2(i+1)-ind2(i)+1;
        if sub<nums(i)
            ind2(i+1) = ind2(i+1)+(nums(i)-sub);
        end
    end


   % Generate blocks amplitude randomly
   rho=[];gen_B=[];gen_G=[];gen_Cov=[];blks={};
  
    for j=1:blkNum
        rho(j) = rand(1,1); 
        gen_B{j} = zeros(nums(j));
        eleB = 1; eleG = 1; 
        for i=1:(nums(j)-1)
            eleB = [eleB , rho(j)^i];
            eleG = [eleG , randi(20)];
        end
        gen_B{j} = toeplitz(eleB);
        gen_G{j} = diag(eleG);
        gen_Cov{j} = gen_G{j}*gen_B{j}*gen_G{j};

        % mean and covariance
        custom_mean = zeros(nums(j),1); 
        custom_covariance = gen_Cov{j}; 
        
        % multi-normal variable vector
        num_samples = 1; 
        random_numbers = mvnrnd(custom_mean, custom_covariance, num_samples);

        blks{j} = random_numbers;
    end
    
    % put blocks at locations
    Wgen = zeros(N,1);  blkLoc2 = [];
    for i = 1:blkNum
        nn = ind2(i)+nums(i)-1;
        mark = length(blkLoc2);
        blkLoc2 = [blkLoc2,ind2(i):(nn)];
        Wgen(ind2(i):(nn)) = blks{i};
    end
    if length(Wgen)>N
        Wgen=Wgen(1:end-(length(Wgen)-N));
    end

%% Noisy Measurements

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
        
    
    
%% ================== BSBL-EM  ==============
    blkStartLoc = [1:h_em_l:N];
    learnlambda = 1;         % when SNR < 25dB, set learnlambda = 1;
                             % when SNR >=25dB, set learnlambda = 2;
                             % when no noise present, set learnlambda = 0;
    tic;
    
    Result1 = BSBL_EM(Phi,Y,blkStartLoc,learnlambda);

    t_em(it) = toc;
    mse_em(it) = (norm(Wgen - Result1.x(:,Result1.count),'fro')/norm(Wgen,'fro'))^2;


    fprintf('BSBL-EM(learn correlation) : time: %4.3f, MSE: %g\n',mean(t_em),mean(mse_em));



%% ================== DivSBL  ==============
    blkStartLoc_l = [1:h_em_l:N];
    tic;
    Result_DV = DivSBL(Phi, Y, blkStartLoc_l);
    
    t_em_DV(it) = toc;
    mse_em_DV(it) = (norm(Wgen - Result_DV.x(:,Result_DV.count),'fro')/norm(Wgen,'fro'))^2;


    fprintf('DivSBL: time: %4.3f, MSE: %g\n',mean(t_em_DV),mean(mse_em_DV));

%% ================== PC-SBL  ==============

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


%% ================== SBL(RVM)  ==============
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

%% ================== Group Lasso  ==============
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

%% ================== Group BPDN  ==============
sigma = 1e-6;
tic;
numElements = N / h_em_l;
repeatedArray = repmat(1:numElements, h_em_l, 1);
groups = repeatedArray(:);
x_GBPDN = spg_group( A, Y, groups, sigma);
t_GBPDN(it) = toc;

mse_GBPDN(it) = (norm(Wgen - x_GBPDN,'fro')/norm(Wgen,'fro'))^2;


fprintf('Group BPDN: time: %4.3f, MSE: %g\n',mean(t_GBPDN),mean(mse_GBPDN));


%% ================== StructOMP  ==============
H=N; W=1; n=H*W; P=1; 
K=blkNum*blkLen; %total non-zero num %Pre-set
q=blkNum;% block Num  
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

%% Reconstruction figure
figure(1);%Original

stem(Wgen, 'LineWidth', 1.7, 'Marker', 'none', 'Color', '#0072BD'); 
xlim([min(Wgen), N]);

title('Original', 'FontSize', 14);


%% DivSBL
figure(2);

stem(Result_DV.x(:,Result_DV.count), 'LineWidth', 1.7, 'Marker', 'none', 'Color', '#0072BD'); % 浣跨RGB棰插艰〃绀鸿?
xlim([min(Result_DV.x(:,Result_DV.count)), N]);

title('DivSBL', 'FontSize', 14); 


%% BSBL-EM
figure(3);

stem(Result1.x(:,Result1.count), 'LineWidth', 1.7, 'Marker', 'none', 'Color', '#0072BD'); % 浣跨RGB棰插艰〃绀鸿?
xlim([min(Result1.x(:,Result1.count)), N]);

title('BSBL-EM', 'FontSize', 14); 



%% %PC-SBL
figure(4);

stem(x_new, 'LineWidth', 1.7, 'Marker', 'none', 'Color', '#0072BD'); % 浣跨RGB棰插艰〃绀鸿?
xlim([min(x_new), N]);

title('PC-SBL', 'FontSize', 14); 

%% %SBL
figure(5);

stem(x_BCS, 'LineWidth', 1.7, 'Marker', 'none', 'Color', '#0072BD'); % 浣跨RGB棰插艰〃绀鸿?
xlim([min(x_BCS), N]);
%% %Group Lasso
figure(6);

stem(x_group_lasso, 'LineWidth', 1.7, 'Marker', 'none', 'Color', '#0072BD'); % 浣跨RGB棰插艰〃绀鸿?
xlim([min(x_group_lasso), N]);

title('Group Lasso', 'FontSize', 14); 

%% % Group BPDN
figure(7);

stem(x_GBPDN, 'LineWidth', 1.7, 'Marker', 'none', 'Color', '#0072BD'); % 浣跨RGB棰插艰〃绀鸿?
xlim([min(x_GBPDN), N]);

title('Group BPDN', 'FontSize', 14); 

%% % StructOMP
figure(8);

stem(x_StrOMP, 'LineWidth', 1.7, 'Marker', 'none', 'Color', '#0072BD'); % 浣跨RGB棰插艰〃绀鸿?
xlim([min(x_StrOMP), N]);

title('StructOMP', 'FontSize', 14); 

%% Reference
%[1]Zhang, Z. and Rao, B. D. Extension of sbl algorithms for the recovery of block sparse signals with intra-block correlation. IEEE Transactions on Signal Processing, 61(8):2009?015, 2013.

