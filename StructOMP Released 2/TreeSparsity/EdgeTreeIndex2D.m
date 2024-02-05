function [Edge2, Edge4]=EdgeTreeIndex2D(H,W, IdxParent, IdxChildren);
% if length(IdxParent)~=H
%     fprintf('error in Dim');
% end
N=length(IdxParent);
[Cn, Cn2]=size(IdxChildren); 
Is=[1:N]';

Isp=IdxParent;

Isc=zeros(N, Cn2); 
Isc(1:Cn,1:Cn2)=IdxChildren(:,:);

Edge2=[Is, Isp];

Edge4=[Is, Isp,Isc];
return