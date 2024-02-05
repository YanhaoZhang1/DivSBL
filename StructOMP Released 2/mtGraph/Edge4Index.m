function Edge4=Edge4Index(H, W);
N = H*W;
Edge4=zeros(N,5);
%% Self
Edge4(:,1)=[1:N;];
%% down
is=Edge4(:,1)+1;
is(H:H:N,1)=[H-1:H:N-1;];
Edge4(:,2)=is;
%% up
is=Edge4(:,1)-1;
is(1:H:N-H+1,1)=[2:H:N-H+2;];
Edge4(:,3)=is;
%% LEFT
is=Edge4(:,1)-H;
is(1:H,1)=[H+1:2*H;];
Edge4(:,4)=is;
%% RIGHT
is=Edge4(:,1)+H;
is(N-H+1:N,1)=[N-2*H+1:N-H;];
Edge4(:,5)=is;
return