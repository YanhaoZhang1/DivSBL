function Result = DivSBL(Phi, y, blkStartLoc, varargin)
%%% Funtion for DivSBL.
%
%   The algorithm solves the inverse problem for the block sparse model:
%                        y = Phi * x + n
%
%
% ============================== INPUTS ============================== 
%   Phi         : N X M known matrix
%
%   y           : N X 1 measurement vector 
%
%   blkStartLoc : Start location of each block
%                 
% [varargin values -- in most cases you can use the default values]
%
%  'PRUNE_GAMMA'  : threshold to prune out small gamma_i 
%                   (generally, 10^{-3} or 10^{-2})
%
%  'LAMBDA'       : user-input value for lambda
%                  [ Default: LAMBDA=1e-14 when LearnLambda=0; LAMBDA=std(y)*1e-2 in noisy cases ]
%
%  'MAX_ITERS'    : Maximum number of iterations.
%                 [ Default value: MAX_ITERS = 800 ]
%
%  'EPSILON'      : Solution accurancy tolerance parameter 
%                 [ Default value: EPSILON = 1e-8   ]
%
%  'PRINT'        : Display flag. If = 1: show output; If = 0: supress output
%                 [ Default value: PRINT = 0        ]
%
% ==============================  OUTPUTS ============================== 
%   Result : A structured data with:
%      Result.x          : the estimated block sparse signal
%      Result.gamma_used : indexes of nonzero groups in the sparse signal
%      Result.gamma_est  : the gamma values of all the groups of the signal
%      Result.B0          : the final value of the diviersifed correlation B
%      Result.count      : iteration times
%      Result.lambda     : the final value of lambda
%
%
%
%    < Full-Command Example >
%           Result =  BSBL_EM(Phi, y, blkStartLoc,                                              
%                                                 'lambda', 1e-3,...
%                                                 'prune_gamma',1e-2,...
%                                                 'MAX_ITERS', 800,...
%                                                 'EPSILON', 1e-8,...
%                                                 'PRINT',0);
%
%
% ============= Version =============
%   1.0 (01/27/2024)
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scaling...
scl = std(y);
if (scl < 0.4) | (scl > 1)
    y = y/scl*0.4;
end


% Default Parameter Values for Any Cases
EPSILON       = 1e-4;       % solution accurancy tolerance  
MAX_ITERS     = 800;        % maximum iterations
PRINT         = 0;          % don't show progress information

lambda = 1e-3;    
% PRUNE_GAMMA = 1*1e-2;
PRUNE_GAMMA = 0.6*1e-2;   %for Images and audio


 

if(mod(length(varargin),2)==1)
    error('Optional parameters should always go by pairs\n');
else
    for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'prune_gamma'
                PRUNE_GAMMA = varargin{i+1}; 
            case 'lambda'
                lambda = varargin{i+1};    
            case 'epsilon'   
                EPSILON = varargin{i+1}; 
            case 'print'    
                PRINT = varargin{i+1}; 
            case 'max_iters'
                MAX_ITERS = varargin{i+1};  
            otherwise
                error(['Unrecognized parameter: ''' varargin{i} '''']);
        end
    end
end



if PRINT
    fprintf('\n====================================================\n');
    fprintf('           Running DivSBL....... \n');
    fprintf('           Information about parameters...\n');
    fprintf('====================================================\n');
    fprintf('PRUNE_GAMMA  : %e\n',PRUNE_GAMMA);
    fprintf('lambda       : %e\n',lambda);   
    fprintf('EPSILON      : %e\n',EPSILON);
    fprintf('MAX_ITERS    : %d\n\n',MAX_ITERS);
end


%% Initialization
[N,M] = size(Phi);
Phi0 = Phi;
blkStartLoc0 = blkStartLoc;
p = length(blkStartLoc);   % block number
for k = 1 : p-1
    blkLenList(k) = blkStartLoc(k+1)-blkStartLoc(k);
end
blkLenList(p) = M - blkStartLoc(end)+1;
maxLen = max(blkLenList);

for k = 1 : p
    Sigma0{k} = eye(blkLenList(k)); %The initial Cov matrix for each block
end

gamma = ones(M,1);
keep_list = [1:p]';
usedNum = length(keep_list);
mu_x = zeros(M,1);
count = 0;


%% Iteration
while (1)
    
    count = count + 1;
    sub_gamma=[];
    currentLoc = 0;

    %=========== Prune weights as their hyperparameters go to zero ==============
    for i = 1:usedNum
        currentLen = size(Sigma0{i},1);
        currentLoc = currentLoc + 1;
        segp{i} = currentLoc : 1 : currentLoc + currentLen - 1;
        sub_gamma(i) = mean(gamma(segp{i}));
        tt = cell2mat(segp(:,i));
        currentLoc = tt(end);
    end

    if (min(sub_gamma) < PRUNE_GAMMA)
        index = find(sub_gamma > PRUNE_GAMMA);
        I = cell2mat(segp(:,index));
        usedNum = length(index);
        keep_list = keep_list(index);
        if isempty(keep_list), 
            fprintf('\n====================================================================================\n');
            fprintf('x becomes zero vector. The solution may be incorrect. \n');
            fprintf('Current ''prune_gamma'' = %g, and Current ''EPSILON'' = %g.\n',PRUNE_GAMMA,EPSILON);
            fprintf('Try smaller values of ''prune_gamma'' and ''EPSILON'' or normalize ''y'' to unit norm.\n');
            fprintf('====================================================================================\n\n');
            break; 
        end;
        blkStartLoc = blkStartLoc(index);
        blkLenList = blkLenList(index);
        
        % prune gamma and associated components in Sigma0 
        sub_gamma = sub_gamma(index);
        gamma=gamma(I);   
        temp = Sigma0;
        Sigma0 = [];
        for k = 1 : usedNum
            Sigma0{k} = temp{index(k)};
        end
        
        % construct new Phi
        temp = [];
        for k = 1 : usedNum
            temp = [temp, Phi0(:,blkStartLoc(k):blkStartLoc(k)+blkLenList(k)-1)];
        end
        Phi = temp;
        %clear temp;
    end

    %=================== Compute new weights =================
    mu_old = mu_x;
    PhiBPhi = zeros(N);
    currentLoc = 0;
    for i = 1 : usedNum
        
        currentLen = size(Sigma0{i},1);
        currentLoc = currentLoc + 1;
        currentSeg = currentLoc : 1 : currentLoc + currentLen - 1;
        
        PhiBPhi = PhiBPhi + Phi(:, currentSeg)*Sigma0{i}*Phi(:, currentSeg)'; 
        currentLoc = currentSeg(end);
    end
    
    H = Phi' /(PhiBPhi + lambda * eye(N));
    Hy = H * y;      
    HPhi = H * Phi;
    
    mu_x = zeros(size(Phi,2),1);
    Sigma_x = [];
    Cov_x = [];
    G = []; GCG=[]; invG=[];
     
    B = []; invB = []; B0 = zeros(maxLen); r0 = zeros(1); r1 = zeros(1);
    currentLoc = 0;
    
    for i = 1 : usedNum
        
        currentLen = size(Sigma0{i},1);
        currentLoc = currentLoc + 1;
        seg = currentLoc : 1 : currentLoc + currentLen - 1;
        
        mu_x(seg) = Sigma0{i} * Hy(seg);       % solution
        Sigma_x{i} = Sigma0{i} - Sigma0{i} * HPhi(seg,seg) * Sigma0{i};
        Cov_x{i} = Sigma_x{i} + mu_x(seg) * mu_x(seg)';
        currentLoc = seg(end);

        G{i} = diag(sqrt(gamma(seg)));
        invG{i}=inv(G{i});
        GCG{i}=invG{i}*Cov_x{i}*invG{i};
        B0 = B0 + GCG{i};
    end
    B0=B0/usedNum;
    
    % Fist Constrain all the blocks have the same correlation structure
    % (an effective strategy to avoid overfitting)
    %%%%%%%%%%  Common B %%%%%%%%%%%%%%%%%%
    b = (mean(diag(B0,1))/mean(diag(B0)));
    if abs(b) >= 0.99, b = 0.98*sign(b); end;
    bs = [];
    for j = 1 : maxLen, bs(j) = (b)^(j-1); end;
    B0 = toeplitz(bs);

    for i = 1 : usedNum

        B{i} = B0;
        invB{i} = inv(B{i});
    end

  %%%%%%%%%%%%  Diversified Bi %%%%%%%%%%%%%%%%%%
    % Dual Acsent
    lam=zeros(usedNum,1);
    eps=1e-3;   
    for i = 1 : usedNum

            alpha=0.01;
            B{i}=(GCG{i})/(1+2*lam(i));
            lam(i)=lam(i)+alpha*(log(det(B{i}))-log(det(B0)));

        %Toepliz correction?            
        b = (mean(diag(B{i},1))/mean(diag(B{i})));
        if abs(b) >= 0.99, 
            b = 0.98*sign(b); 
        end;
        bs = [];
        for j = 1 : maxLen, bs(j) = (b)^(j-1); end;
        B{i} = toeplitz(bs);

        invB{i} = inv(B{i});
    end

%%%%%%%%%%%% Estimate gamma(i) and lambda  %%%%%%%%%%%% 
    gamma_old = gamma;
    lambdaComp = 0; currentLoc = 0;
    for i =  1 : usedNum

        currentLen = size(Sigma_x{i},1);
        currentLoc = currentLoc + 1;
        currentSeg = currentLoc : 1 : currentLoc + currentLen - 1;
        lambdaComp = lambdaComp + trace(Phi(:,currentSeg)*Sigma_x{i}*Phi(:,currentSeg)');
        currentLoc = currentSeg(end);

        DiaginvB = diag(invB{i});             
        DiagCov = diag(Cov_x{i});
        A = DiaginvB.*DiagCov;

        for j = 1:currentLen
            bb = invB{i}(j,:);                
            gg = diag(1./G{i})';
            gg(j) = 0; 
            bg = bb.*gg;
            T = bg*Cov_x{i}(:,j);
            xx1 = (-T+sqrt(T^2+4*A(j)))/(2*A(j));
            G{i}(j,j) = 1/xx1; 
        end
        gamma(currentSeg) = (diag(G{i})).^2;
        Sigma0{i} = G{i}*B{i}*G{i};
    end
    lambda = norm(y - Phi * mu_x,2)^2/N + lambdaComp/N;    
 
    
    % ================= Check stopping conditions, eyc. ==============
    if (size(mu_x) == size(mu_old))
        dmu = max(max(abs(mu_old - mu_x)));
        if (dmu < EPSILON)  
            break; 
        end;
    end;
    if (PRINT) 
        disp([' iters: ',num2str(count),...
            ' num coeffs: ',num2str(usedNum), ...
            ' min gamma: ', num2str(min(gamma)),...
            ' gamma change: ',num2str(max(abs(gamma - gamma_old))),...
            ' mu change: ', num2str(dmu)]); 
    end;
    if (count >= MAX_ITERS), 
%         fprintf('count >= MAX_ITERS');if PRINT, fprintf('Reach max iterations. Stop\n\n'); end; 
        break;  
    end;
%% Expand hyperparameyers
sub_gamma_used = sort(keep_list);
I = cell2mat(segp(:,sub_gamma_used));
gamma_est(I) = gamma;
lam_est = zeros(p,1);

%% reconstruct the original signal
x1 = zeros(M,1);
currentLoc = 0;
for i = 1 : usedNum

    currentLen = size(Sigma0{i},1);
    currentLoc = currentLoc + 1;
    seg = currentLoc : 1 : currentLoc + currentLen - 1;

    realLocs = blkStartLoc0(keep_list(i)) : blkStartLoc0(keep_list(i))+currentLen-1;

    x1( realLocs ) = mu_x( seg );
    currentLoc = seg(end);
end

if (scl < 0.4) | (scl > 1)
    Result.x(:,count) = x1 * scl/0.4;
else
    Result.x(:,count) = x1;
end

Result.B = B;
Result.count = count;
Result.lambda = lambda;
Result.multi(:,count) = lam_est;
b = (mean(diag(B0,1))/mean(diag(B0)));
if abs(b) >= 0.99, b = 0.98*sign(b); end;
bs = [];
for j = 1 : maxLen, bs(j) = (b)^(j-1); end;
B0 = toeplitz(bs);
Result.B0 = B0;
Result.gamma(:,count) = gamma_est;
end;

return;

