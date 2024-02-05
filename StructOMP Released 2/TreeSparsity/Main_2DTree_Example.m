%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This script is used to obtain results for strong line structured sparsity 
%%%% Junzhou Huang, Tong Zhang, Dimitris Metaxas "Learning with Structured Sparsity"
%%%% By Junzhou Huang, Rutgers University,jzhuang@cs.rutgers.edu
%%%% Jan., 2009    Updated  Dec. 20, 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all
load ImgData;
x=ImgData;
rng(1,'twister')
[H, W, P]=size(x);
figure; imshow(x, []);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D wavelet transformation and sparse signal -- theta0 
DecLevel=6;
WaveletName='db1';
[C0, S0] = wavedec2(x, DecLevel, WaveletName);
scaling0=C0(1:S0(1,1)*S0(1,2));
theta0=C0(S0(1,1)*S0(1,2)+1:end)';

[IdxParent, IdxChildren, Ms]=WaveRelation2D_NoSc(C0, S0);
[H, W]=size(theta0);
plot(theta0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=length(theta0);
[Edge2,Edge4]=EdgeTreeIndex2D(H, W, IdxParent', IdxChildren');
[B, Bm]=GetTreeBlocksMatrix(H, W, IdxParent,2);
[BC, BCm]=GetBlocksConnectionMatrix(Bm, Edge2);   
save BCm BCm; save Bm Bm;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Bm; load BCm;


[sv1,si1]=sort(abs(theta0));
sv2=cumsum(sv1)/sum(sv1);
sind2=find(sv2<=0.05);%%%%%%%%%%%
index=si1(sind2);
index2=TraceChildNode(index, Edge4);

theta=theta0; 
% theta(si1(sind3))=0;
theta(index2)=0;
K=length(find(theta~=0));

mq=length(Bm); q=sum(theta(1:S0(2,1)*S0(2,2))~=0);
p=length(theta0);

m_num=2048;

for j=1:length(m_num),
    fprintf('j=%d\n', j);
    for iter=1:1
        p=length(theta0);
        m=round(m_num(j));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% step 2: Creating random projection matrix and measurements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=1:P
            A=randn(m,p);%/A=single(A);
            item=sqrt(sum(A.*A,2));
            A=A./repmat(item, [1,p]);
            AA{i}=A; e=randn(m,1); 
            y0=AA{i}*theta0;
            norm_y=1;
            y{i}=y0+20*e;
        end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% step 3: Sparse recovery
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
       %%% method 1: OMP
        tic;
%         [ta1, Rec1] = OMP_CS(K,AA{1},y{1}, 1);
        [ta1, Rec1] = OMPlscx0_CS(K,AA{1},y{1});
        ta1=ta1*norm_y;
        time1=toc;
        xest1=waverec2([scaling0,ta1'], S0, WaveletName);
        fprintf('method OMP\n');
        
        %%% method 2: Lasso
        tt=cputime;
        [ta2, numIters, activationHist, duals] = SolveLasso(AA{1}, y{1}, p, 'lasso', 1200);    
        ta2=ta2*norm_y;
        time2=cputime-tt;
        xest2=waverec2([scaling0,ta2'], S0, WaveletName);
        fprintf('method Lasso\n')
    
        %%% method 2: StructOMP
        tic;
        lamada=1;
        cl0=lamada*q*log2(mq)+K;
        [ta3, input, norm_save3] = GraphOMP_CS(cl0,AA{1},y{1}, Bm, BCm, lamada, []);
        ta3=ta3*norm_y;
        time3=toc;
        xest3=waverec2([scaling0,ta3'], S0, WaveletName);
        fprintf('method TreeStructOMP\n');
        
        diff1=norm(xest1(:)-x(:),2)/norm(x(:),2);
        diff2=norm(xest2(:)-x(:),2)/norm(x(:),2);  
        diff3=norm(xest3(:)-x(:),2)/norm(x(:),2);  

        result(:,iter)=[diff1,diff2,diff3, time1 time2, time3];
        figure;
        subplot(1,4,1); imshow(x, []); title('(a)');
        subplot(1,4,2); imshow(xest1, []); title('(b)');
        subplot(1,4,3); imshow(xest2, []); title('(c)');
        subplot(1,4,4); imshow(xest3, []); title('(d)')
        
        %%% method 3: BSBL
        
            %================== BSBL-EM exploiting intra-block correlation ==============
    h_em_l = 45;
    Phi=AA{1}; Y=y{1};
    blkStartLoc = [1:h_em_l:N];
%     blkStartLoc = [1:blkLen:N];
    learnlambda = 1;         % when SNR < 25dB, set learnlambda = 1;
                             % when SNR >=25dB, set learnlambda = 2;
                             % when no noise present, set learnlambda = 0;
    tic;
    
    Result1 = BSBL_EM(Phi,Y,blkStartLoc,learnlambda);
    t_em = toc;
    xest_BSBL= waverec2([scaling0,Result1.x(:,Result1.count)'], S0, WaveletName);

    mse_em = norm(xest_BSBL(:)-x(:),2)/norm(x(:),2);

    fprintf('BSBL-EM(learn correlation) : time: %4.3f, MSE: %g\n',mean(t_em),mean(mse_em));

%       =========================================================
%                    Diver-Gamma-BSBL-EM with h = 6
%     =========================================================

    blkStartLoc_l = [1:h_em_l:N];
    learnlambda_em_l = 1;          % when SNR is low (e.g. < 20dB), set to 1
%                                   when SNR is high (e.g. <20dB), set to 2
%                                   in noiseless cases, set to 0
    tic;
    Result_DV = Div_Var_BSBL_EM(Phi, Y, blkStartLoc_l, learnlambda_em_l);
    t_em_DV= toc;
   xest_DivSBL= waverec2([scaling0,Result_DV.x(:,Result_DV.count)'], S0, WaveletName);

    mse_DV = norm(xest_DivSBL(:)-x(:),2)/norm(x(:),2);
    fprintf('DivSBL : time: %4.3f, MSE: %g\n',mean(t_em_DV),mean(mse_DV));

end
CompareResults(:,j)=mean(result, 2)
CompareResultsStd(:,j)=std(result,0, 2);
end