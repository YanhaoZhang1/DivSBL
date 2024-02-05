function Rec = Model_CS(K, A, y, Jb, beta, H, W, Option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Model based (Compressive Sensing) for overlaped group sparsity
%%% Input:
%%%   K:   the sparsity  group number
%%%   A:   projection matrix. Cell structure, One cell for one task.
%%%   y:   CS measurements. Cell structure, One cell for one task.
%%%   Jb:  the block size for each task
%%%   beta: weight matrix for neighbors of each sparse coefficient
%%%   [H, W]:  for 1D data, W=1; for 2D data, they are the sizes
%%%   Option: 'None' not group clustering prior, beta=0
%%%           'Grid' 4 neighbors in 2D grid
%%%           'GridCol'  2 neighbors in each column. No relation within rows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Variable:
%%%  Jy =length(y), the number of tasks. 
%%%  [dJy, TJb]=size(A{Jy})
%%%  Jx=Jb*Jy, the total group size for all tasks
%%%  Jx*K: Total sparsity number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Output:
%%%   Rec.x_hat [TJb, Jy]: the reconstructed sparse variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Rutgers University, Nov 6, 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iscell(y)
    Jy = length(y);
else
    Jy = 1; A = {A}; y = {y};
end
Jx=Jy*Jb; TJb=size(A{1}, 2); 
T=TJb/Jb; TJx=T*Jx;
y_r = y; in = 1;
for jj=1:Jy,
    cv(1, (jj-1)*TJb+1:jj*TJb) = y_r{jj}'*A{jj};
end
index_T=[1:T]';
cv_index=S_Prune(cv', index_T, Jx, Jb, K, beta, H, W, Option);
Index_save(:,in) = cv_index;
for jj=1:Jy,
    A_x{jj} = A{jj}(:, cv_index);
    x_p((jj-1)*Jb*K+1:jj*Jb*K, 1) = inv(A_x{jj}'*A_x{jj})*A_x{jj}' * y{jj};
    y_r{jj} = y{jj} - A_x{jj}*x_p( (jj-1)*Jb*K+1:jj*Jb*K, 1);
    norm_yr(jj) =norm(y_r{jj});
end
norm_save(in) = norm(norm_yr);
x_curr=x_p;
while 1
    in = in+1;
    for jj=1:Jy,
        cv(1, (jj-1)*TJb+1:jj*TJb) =  y_r{jj}'*A{jj};
    end

    index_T=[1:T]';
    cv_index=S_Prune(cv', index_T, Jx, Jb, 1*K, beta, H, W, Option);
    a1=reshape(Index_save(:,in-1), [K, Jb]);
    a2=reshape(cv_index, [1*K, Jb]);
    cv_add=unique([a1;a2], 'rows');
    cv_add=cv_add(:)';
    cv_length=length(cv_add);
    x_p=zeros(Jy*cv_length,1);
    for jj=1:Jy,
        A_x{jj} = A{jj}(:, cv_add);
        x_p((jj-1)*cv_length+1:jj*cv_length, 1) = pinv(A_x{jj}'*A_x{jj})*A_x{jj}' * y{jj};
    end   
 
    xpindex_T=cv_add(1,1:cv_length/Jb)';
    i_sort=S_Prune(x_p, xpindex_T, Jx, Jb, K, beta(xpindex_T), H, W, Option);
    cv_index = cv_add( i_sort);
    cv_index = sort( cv_index );
    Index_save(:,in)=cv_index;
    x_p=zeros(Jy*Jb*K,1);
    for jj=1:Jy,
        A_x{jj} = A{jj}(:, cv_index);
        x_p((jj-1)*Jb*K+1:jj*Jb*K, 1) = pinv(A_x{jj}'*A_x{jj})*A_x{jj}' * y{jj};
        y_r{jj} = y{jj} - A_x{jj}*x_p( (jj-1)*Jb*K+1:jj*Jb*K, 1);
        norm_yr(jj) =norm(y_r{jj});
    end
    norm_save(in) = norm(norm_yr);
    if ( norm_save(in) == 0 ) ||(norm_save(in)>=norm_save(in-1) )
        break;
    end
    x_curr=x_p;
end
x_hat = zeros(TJx,1);
sup=[];
for jj=1:Jy,
    temp=Index_save(:,in-1)+TJb*(jj-1);
%     A_x{jj} = A{jj}(:, temp);
%     x_curr_Jy((jj-1)*Jy*Jb*K+1:jj*Jy*Jb*K, 1) = inv(A_x{jj}'*A_x{jj})*A_x{jj}' * y{jj};
    sup=[sup; temp];
end

x_hat(sup,1) = reshape(x_curr,Jx*K,1);
Rec.T = Index_save;
Rec.x_hat = x_hat;
Rec.PResidue = norm_save;
return


function gs_index=S_Prune(x, xindex, Jx, Jb, K, beta, H, W, Option);
%%%% x=[X1; X2; ...; XJ], XJ-> T x 1   
%%%% x [TJx, 1]  xindex [T, 1]
%%%% K is sparse for T, thus total sparsity is K*Jx for T*Jx
%%%% beta is [T x 1]
%%%   gs_index [1 x K*Jb]
TJx=length(x); T=TJx/Jx; x0=x; 

if strcmp(Option, 'None')
        T=TJx/Jx;
        if T~=round(T) disp('GroupSize Error');   end
        XM=reshape(x, [T, Jx]);
        XMJnorm=sqrt(sum(XM.*XM, 2));
        [xmj_sort, xmj_index] = sort(XMJnorm,'descend');
        xmj_index = sort( xmj_index(1:K) );
        gs_index=[];
        for jj=1:Jb,
            temp=xmj_index+T*(jj-1);
            gs_index=[gs_index; temp];
        end
    
elseif strcmp(Option,'Grid')
        Edge4=Edge4Index(H, W);  
        m=max(H*W, T);
        betaM=repmat(beta, [1, Jx]);              
        XM=reshape(x, [T, Jx]);
        xtmp=zeros(m,Jx);
        xtmp(xindex,:)=XM;
        for enum=2:size(Edge4,2)
            tmp=Edge4(xindex,enum);      
            XM=[XM, betaM.*xtmp(tmp,:)];
        end
        XMJnorm=sqrt(sum(XM.*XM, 2));
        [xmj_sort, xmj_index] = sort(XMJnorm,'descend');
        xmj_index = sort( xmj_index(1:K) );
        gs_index=[];
        for jj=1:Jb,
            temp=xmj_index+T*(jj-1);
            gs_index=[gs_index; temp];
        end
elseif strcmp(Option,'GridCol') %% FOR EXAMPLS: H=nc; W=classNum, thus T=nc*classNum;
        Edge4=Edge4Index(H, W);
        Edge2=Edge4(:,1:3);
        m=max(H*W, T);
        betaM=repmat(beta, [1, Jx]);  
        XM=reshape(x, [T, Jx]);
        xtmp=zeros(m,Jx);
        xtmp(xindex,:)=XM;
        for enum=2:size(Edge2,2)
            tmp=Edge2(xindex,enum);      
            XM=[XM, betaM.*xtmp(tmp,:)];
        end
        XMJnorm=sqrt(sum(XM.*XM, 2));
        [xmj_sort, xmj_index] = sort(XMJnorm,'descend');
        xmj_index = sort( xmj_index(1:K) );
        gs_index=[];
        for jj=1:Jb,
            temp=xmj_index+T*(jj-1);
            gs_index=[gs_index; temp];
        end    
end     
return

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