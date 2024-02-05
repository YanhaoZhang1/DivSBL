%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This script is used to obtain results for graph sparsity 
%%%% Junzhou Huang, Tong Zhang, Dimitris Metaxas "Learning with Structured Sparsity"
%%%% By Junzhou Huang, Rutgers University,jzhuang@cs.rutgers.edu
%%%% Jan., 2009    Updated  Dec. 20, 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear; close all;
%%% parameters
load ImgMask1;
[T, WI]=size(ImgMask1);
[H, W, P]=size(ImgMask1{1,1});
n=H*W;
Edge4=Edge4Index(H, W);
[B, Bm]=GetBlocksMatrix(H,W,1);
[BC, BCm]=GetBlocksConnectionMatrix(Bm, Edge4);
mq=size(Bm, 1); clear B; clear BC;

save BCm BCm; save Bm Bm;

% load BCm; load Bm;

m_num=[2 2.5 3 3.5 4];

for j=1:length(m_num),
      fprintf('j=%d\n', j);
    for bgb=1:5,
        fprintf('bgb=%d\n', bgb);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% step 1: Get DG(K, q)-sparse 2D background substracted color images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bg=107;
    mask=ImgMask1{bg, 2}; bin=repmat(mask, [1, 1, P]);
    Img0=ImgMask1{bg, 1};
    x=double(Img0).*bin;
%     figure; imshow(uint8(x), []); 
    
    K=length(find(x(:)>0))/P; 
    m=round(m_num(j)*K); q=4;
%     m=900;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% step 2: Creating random projection matrix and measurements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x=reshape(x, [n, P]);
    for i=1:P
        A=randn(m,n);item=sqrt(sum(A.*A,2));A=A./repmat(item, [1,n]);
        AA{i}=A; e=randn(m,1); %e=e/norm(e(:));        
        y{i}=AA{i}*x(:,i)+0.01*e;
    end
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% step 3: Sparse recovery
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  %%% method 1: OMP    
%    tic;
%     for i=1:P,
% %         [xest, Rec1] = OMP_CS(K,AA{i},y{i}, 1);
%         [xest, Rec1] = OMPlscx0_CS(K,AA{i},y{i});
%         xest1(:,:,i)=reshape(xest, [H, W]);
%     end
%     time1=toc;
    %%% method 2: Lasso
    tic;
    for i=1:P,
        [xest, numIters] = SolveLasso(AA{i}, y{i}, n, 'lasso', 1000);
        xest2(:,:,i)=reshape(xest, [H, W]);
    end    
    time2=toc;  
%    
%    %%% method 1: StructOMP
%     tic;
%     for i=1:P,
%        lamada=1; cl0=lamada*q*log2(mq)+K;
%        [xest, input1, norm_save1] = GraphOMP_CS(cl0,AA{i},y{i}, Bm, BCm, lamada, []);
%        xest3(:,:,i)=reshape(xest, [H, W]);
%     end    
%     time3=toc; 
%     %%%% Normalized Recovery Error
%     diff1=norm(xest1(:)-x(:),2)/norm(x(:),2);
    diff2=norm(xest2(:)-x(:),2)/norm(x(:),2);
%     diff3=norm(xest3(:)-x(:),2)/norm(x(:),2);
   
    result(:,bgb)=[diff2,diff2 diff2, time2, time2  time2];
%     result(:,bgb)=[diff1, diff3, time1,  time3];
end
CompareResults(:,j)=mean(result, 2)
CompareResultsStd(:,j)=std(result,0, 2);
end 
CompareResults2=CompareResults;
CompareResultsStd2=CompareResultsStd;
save CompareResults2 CompareResults2;
save CompareResultsStd2 CompareResultsStd2;

figure; hold on;
errorbar(m_num*K, CompareResults(1,:), CompareResultsStd(1,:), 'bv-', 'linewidth', 2); 
errorbar(m_num*K, CompareResults(2,:), CompareResultsStd(2,:), 'gd-', 'linewidth', 2); 
errorbar(m_num*K, CompareResults(3,:), CompareResultsStd(3,:), 'rp-', 'linewidth', 2); 
ylabel('Recovery Error'); xlabel('Sample Size');box on;
legend('OMP', 'Lasso','StructOMP');
% axis([m_num(1)/K-0.5 m_num(end)/K+0.5 -0.1 1.6])

figure; hold on;
errorbar(m_num*K, CompareResults(4,:), CompareResultsStd(4,:), 'bv-', 'linewidth', 2); 
errorbar(m_num*K, CompareResults(5,:), CompareResultsStd(5,:), 'gd-', 'linewidth', 2); 
errorbar(m_num*K, CompareResults(6,:), CompareResultsStd(6,:), 'rp-', 'linewidth', 2); 
ylabel('Time'); xlabel('Sample Size');box on;
legend('OMP', 'Lasso','StructOMP');
% axis([m_num(1)/K-0.5 m_num(end)/K+0.5 -0.1 1.6])