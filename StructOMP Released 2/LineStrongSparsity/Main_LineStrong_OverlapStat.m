%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This script is used to obtain results for strong line structured sparsity
%%%% Compared with OverlapLasso and Model-CS
%%%% Junzhou Huang, Tong Zhang, Dimitris Metaxas "Learning with Structured Sparsity"
%%%% By Junzhou Huang, Rutgers University,jzhuang@cs.rutgers.edu
%%%% Jan., 2009    Updated  Dec. 20, 2009
%%%% MATLAB 2008a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear; close all;
H=512; W=1; n=H*W; P=1; K=64; q=4; 
Edge4=Edge4Index(H, W);
if W==1
    Edge4=Edge4(:,1:3);
end

[B, Bm]=GetBlocksMatrix(H, W, 1);
[BC, BCm]=GetBlocksConnectionMatrix(Bm, Edge4);

mq=length(B);
% m_num=round(K*[2 2.5 3 3.5]);
m_num=round(K*[2 2.5 3 3.5 4 4.5 5]);
for j=1:length(m_num),
    m=m_num(j)
for iter=1:20, %%% adjust, how many time do you want to run experiments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% step 1: randomly create strong line sparse 1D-data
    %%% randomly distribute the K nozero items into q groups
    %%% randomly choose location of each group and keep them disjoined
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x=zeros(n,1);[H, W]=size(x);
    r=abs(randn(q,1)); r=r+1; r=round(r*K/sum(r)); r(q)=K-sum(r(1:q-1));
    g=round(r*n/K); g(q)=n-sum(g(1:q-1));
    g_cum=cumsum(g);
    for i=1:q,        
        loc=randperm(g(i)-r(i)-2);        
        x_tmp=zeros(g(i), 1);
        x_tmp(loc(1)+2:loc(1)+1+r(i))=sign(randn(r(i), 1));
        x(g_cum(i)-g(i)+1:g_cum(i), 1)=x_tmp;        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% step 2: Creating random projection matrix and measurements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:P
        A=randn(m,n);item=sqrt(sum(A.*A,2));A=A./repmat(item, [1,n]);
        AA{i}=A; e=randn(m,1); 
        y{i}=AA{i}*x(:)+0.01*e;
    end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% step 3: Sparse recovery
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      %% method 1: OverlapLasso
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic; 
    tau=0.01*max(abs(AA{1}'*y{1}));tolA = 1.e-2;
    varGroupMatrix=CreatGroupMatrix(H, W, 1);
    [xest] = OverlapLasso(AA{1},y{1},varGroupMatrix,tau);
    xest1   = xest(1:n,1);
    time1=toc;   
  
    %%%%  Method 2 : Model-CS
    tic;
    Rec = Model_CS(K, AA, y, 1, 0.5*ones(n,1), H, W, 'GridCol');
    time2=toc; xest=Rec.x_hat; xest2=xest(1:n,1);
    
    %%%%  Method 3 : StructOMP
     tic;
    lamada=1;
    cl0=2*lamada*q*log2(mq)+K;
    [xest3, input, norm_save3] = GraphOMP_CS(cl0,AA{1},y{1}, Bm, BCm, lamada, []);
    time3=toc;
 
       
    diff1=norm(xest1(:)-x(:),2)/norm(x(:),2);    
    diff2=norm(xest2(:)-x(:),2)/norm(x(:),2);  
    diff3=norm(xest3(:)-x(:),2)/norm(x(:),2);  
%     diff1=diff2; time1=time2;
    
    result(:,iter)=[diff1,diff2,diff3,time1 time2, time3];
end
CompareResults(:,j)=mean(result, 2);
CompareResultsStd(:,j)=std(result,0, 2);
end

% save CompareResults CompareResults;
% save CompareResultsStd CompareResultsStd;

ms=5; ts=12;

figure; hold on;
errorbar(m_num/K, CompareResults(1,:), CompareResultsStd(1,:), 'bv-', 'linewidth', 2,  'MarkerSize',ms); 
errorbar(m_num/K, CompareResults(2,:), CompareResultsStd(2,:), 'gd-', 'linewidth', 2,  'MarkerSize',ms); 
errorbar(m_num/K, CompareResults(3,:), CompareResultsStd(3,:), 'rp-', 'linewidth', 2,  'MarkerSize',ms); 
ylabel('Recovery Error'); xlabel('Sample Size Ratio (n / k)');box on;
legend('OverLasso','ModelCS','StructOMP');
axis([m_num(1)/K-0.5 m_num(end)/K+0.5 -0.1 1.6])
textobj = findobj('type', 'text');
set(textobj, 'fontsize', ts);
h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',ts); 
h_xlabel = get(gca,'YLabel');
set(h_xlabel,'FontSize',ts); 



figure; hold on;
errorbar(m_num/K, CompareResults(4,:), CompareResultsStd(4,:), 'bv-', 'linewidth', 2,  'MarkerSize',ms); 
errorbar(m_num/K, CompareResults(5,:), CompareResultsStd(5,:), 'gd-', 'linewidth', 2,  'MarkerSize',ms); 
errorbar(m_num/K, CompareResults(6,:), CompareResultsStd(6,:), 'rp-', 'linewidth', 2,  'MarkerSize',ms); 
ylabel('Time'); xlabel('Sample Size Ratio (n / k)');box on;
legend('OverLasso','ModelCS','StructOMP');
textobj = findobj('type', 'text');
set(textobj, 'fontsize', ts);
h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',ts); 
h_xlabel = get(gca,'YLabel');
set(h_xlabel,'FontSize',ts); 


