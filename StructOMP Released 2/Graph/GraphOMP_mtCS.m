function [xest, input, norm_save] = GraphOMP_mtCS(cl0,Phi,y, Bm, BC, lamada, input);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Multi-task Graph structured sparsity
%%% Input:
%%%   cl0:  the coding complexity
%%%   Phi:  projection matrix
%%%   y:   measurements  
%%%   Bm:  the defined block matrix, where the i-th row denotes all index
%%%        entires in sparse x included in the i-th block.
%%%   BC:  the connected relations between blocks, where the i-th row
%%%        denotes all index of blocks connected to the i-th blocks 
%%%   lambda: weights in the coding complexity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Output:
%%%   xest: the unknow structured sparse data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   Junzhou Huang, Tong Zhang, Dimitris Metaxas "Learning with
%%%%   Structured Sparsity", Rutgers University.
%%%%   By Junzhou Huang, jzhuang@cs.rutgers.edu
%%%%   Jan 2009, Updated Dec 20, 2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mt=length(y);
[n,p]=size(Phi{1});
if isempty(input)    
    input.supp = [];input.Blist=[];
    input.cl=0; input.RegionNum=0; 
    input.RegionIndex=[];
    for tt=1:mt
        y_r(:,tt) = y{tt};
        x_p{tt}=[];
    end    
else
    for tt=1:mt,
        y_r(:,tt)=y{tt} - Phi{tt}*input.x(:,tt);
        x_p{tt}=input.x_p{tt};
    end
    
end

OptTol = 1e-8;maxiter_lsqr=3;

in = 0;
while in<10000
   
   in = in+1;
   for tt=1:mt,
       cv1(tt,:) = abs( y_r(:,tt)'*Phi{tt} );
   end
   cv=sqrt(sum(cv1.*cv1,1));
%    [input, newid]=PruneBlocktmp(cv, Bm, BC, input, lamada);
   [input.supp, input.Blist, input.cl, input.RegionNum,input.RegionIndex, newid] = GraphPruneMex(cv, Bm, BC, lamada,input.supp, input.Blist, input.cl, input.RegionNum, input.RegionIndex);
   curr_index=input.supp;
   
%    x0=zeros(length(curr_index),3);
%    idd=~ismembc(curr_index, newid);
%    indx0=find(idd==1);
   for tt=1:mt,
       Phi_x = Phi{tt}(:,curr_index);
       x0=zeros(length(curr_index),1);
       idd=~ismembc(curr_index, newid);
       indx0=find(idd==1);
       x0(indx0)=x_p{tt};
       [x_p{tt}, flag] = lsqr_gp(Phi_x, y{tt}, curr_index', OptTol, maxiter_lsqr, [], [], x0);   
%      x_p = inv(Phi_x'*Phi_x)*Phi_x' * y;
       y_r(:,tt) = y{tt} - Phi_x*x_p{tt};       
   end
   norm_save(in) = norm(y_r(:));
   if input.cl>cl0
       break;
   end
end

xest=zeros(p, mt);
for tt=1:mt
    xest(curr_index,tt)=x_p{tt};
end
input.x=xest;
input.x_p=x_p;
return