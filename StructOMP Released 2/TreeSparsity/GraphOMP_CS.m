function [xest, input, norm_save] = GraphOMP_CS(cl0,Phi,y, Bm, BC, lamada, input);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Graph structured sparsity
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
[n,p]=size(Phi);
if isempty(input)
    y_r = y;
    input.supp = [];input.Blist=[];
    input.cl=0; input.RegionNum=0; 
    input.RegionIndex=[];
    x_p=[];
else
    y_r=y - Phi*input.x;
    x_p=input.x_p;
end

OptTol = 1e-8;maxiter_lsqr=3;

in = 0;
while in<10000
   in = in+1;
   cv = abs( y_r'*Phi );
%    [input, newid]=PruneBlocktmp(cv, Bm, BC, input, lamada);
   [input.supp, input.Blist, input.cl, input.RegionNum,input.RegionIndex, newid] = GraphPruneMex(cv, Bm, BC, lamada,input.supp, input.Blist, input.cl, input.RegionNum, input.RegionIndex);
   curr_index=input.supp;
   Phi_x = Phi(:,curr_index);
   x0=zeros(length(curr_index),1);
   idd=~ismembc(curr_index, newid);
   indx0=find(idd==1);
   x0(indx0')=x_p;   
   [x_p, flag] = lsqr_gp(Phi_x, y, curr_index', OptTol, maxiter_lsqr, [], [], x0);   
%    x_p = inv(Phi_x'*Phi_x)*Phi_x' * y;
   y_r = y - Phi_x*x_p;
   norm_save(in) = norm(y_r);
   if input.cl>cl0
       break;
   end
end
xest=zeros(p, 1);
xest(curr_index,1)=x_p;
input.x=xest;
input.x_p=x_p;
return