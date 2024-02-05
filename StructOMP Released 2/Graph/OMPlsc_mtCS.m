function [xest, norm_save] = OMPlsc_mtCS(K,Phi,y)
mt=length(y);
[m,N]=size(Phi{1});
for tt=1:mt,
    y_r(:,tt) = y{tt};
    x_p{tt}=[];
end
in = 0;
curr_index = [];

OptTol = 1e-8;maxiter_lsqr=2;
Ifull = 1:N;    

while in<K
   in = in+1;
   for tt=1:mt,
       cv1(tt,:) = abs( y_r(:,tt)'*Phi{tt} );
   end
   cv=sqrt(sum(cv1.*cv1,1));
   [iv, id] = max(cv);
   curr_index =unique([curr_index, id]);
   for tt=1:mt
       Phi_x = Phi{tt}(:,curr_index);
       x0=zeros(length(curr_index),1);
       indx0=find(curr_index~=id);
       x0(indx0)=x_p{tt};
       [x_p{tt}, flag] = lsqr_gp(Phi_x, y{tt}, curr_index', OptTol, maxiter_lsqr, [], [], x0);  
%    x_p = inv(Phi_x'*Phi_x)*Phi_x' * y;
       y_r(:,tt) = y{tt} - Phi_x*x_p{tt};
   end
   norm_save(in) = norm(y_r(:));
end
xest=zeros(N, mt);
for tt=1:mt
    xest(curr_index,tt)=x_p{tt};
end
return