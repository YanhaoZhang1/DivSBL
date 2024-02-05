%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This script is used to obtain results for strong line structured sparsity 
%%%% Junzhou Huang, Tong Zhang, Dimitris Metaxas "Learning with Structured Sparsity"
%%%% By Junzhou Huang, Rutgers University,jzhuang@cs.rutgers.edu
%%%% Jan., 2009    Updated  Dec. 20, 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all
load ImgData;
x=ImgData;

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

m_num=[1.25 1.5 1.75 2 2.25 2.5 2.75]*K; 
m_num=round(m_num);

for j=1:length(m_num),
    fprintf('j=%d\n', j);
    for iter=1:3,
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
            y{i}=y0+0.00*e;
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

end
CompareResults(:,j)=mean(result, 2)
CompareResultsStd(:,j)=std(result,0, 2);
end
save CompareResults CompareResults;
save CompareResultsStd CompareResultsStd;



figure; hold on;
errorbar(m_num, CompareResults(1,:), CompareResultsStd(1,:), 'bv-', 'linewidth', 2); 
errorbar(m_num, CompareResults(2,:), CompareResultsStd(2,:), 'gd-', 'linewidth', 2); 
errorbar(m_num, CompareResults(3,:), CompareResultsStd(3,:), 'rp-', 'linewidth', 2); 
ylabel('Recovery Error'); xlabel('Sample Size');box on;
legend('OMP', 'Lasso','StructOMP');
% axis([m_num(1)-round(0.25*K) m_num(end)+round(0.25*K) -0 0.3])

figure; hold on;
errorbar(m_num, CompareResults(4,:), CompareResultsStd(4,:), 'bv-.', 'linewidth', 2);
errorbar(m_num, CompareResults(5,:), CompareResultsStd(5,:), 'gd-', 'linewidth', 2);
errorbar(m_num, CompareResults(6,:), CompareResultsStd(6,:), 'rp-', 'linewidth', 2);
ylabel('Running Time (Second)'); xlabel('Sample Size');box on;
legend('OMP','Lasso', 'StructOMP');
% axis([2 4.5 0 0.1])