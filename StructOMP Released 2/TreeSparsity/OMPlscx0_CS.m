function [xest, norm_save] = OMPlscx0_CS(K,Phi,y)

[m,N]=size(Phi);
y_r = y;in = 0;
curr_index = [];

OptTol = 1e-8;maxiter_lsqr=3;
Ifull = 1:N;    
x_p=[];
while in<K
   in = in+1;
   cv = abs( y_r'*Phi );
   [iv, id] = max(cv);
   curr_index =unique([curr_index, id]);
   Phi_x = Phi(:,curr_index);
   x0=zeros(length(curr_index),1);
   indx0=find(curr_index~=id);
   x0(indx0)=x_p;
   [x_p, flag] = lsqr_gp(Phi_x, y, curr_index', OptTol, maxiter_lsqr, [], [], x0);  
%    x_p = inv(Phi_x'*Phi_x)*Phi_x' * y;
   y_r = y - Phi_x*x_p;
   norm_save(in) = norm(y_r);
end
xest=zeros(N, 1);
xest(curr_index,1)=x_p;
return