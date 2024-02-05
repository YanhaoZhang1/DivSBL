function [B, Bmatrix]=GetBlocksMatrix(H, W, Step);
if W==1  %%% 1D signal
    for i=1:H,
        if i<=Step
            B{i}=[1:i+Step];
        elseif i>=H-Step+1
            B{i}=[i-Step:H];
        else
            B{i}=[i-Step:i+Step];
        end        
    end
    
    
else     %%% 2D signal
    Edge4=Edge4Index(H, W);
    for i=1:length(Edge4)
        B{i}=unique(Edge4(i,:));
    end
end

H=length(B);
maxnum=0;
for i=1:H,
    tmp=length(B{i});
    maxnum=max(maxnum, tmp);
end
Bmatrix=zeros(H, maxnum);
for i=1:H,
    tmp=B{i};
    Bmatrix(i,1:length(tmp))=tmp;
    clear tmp;
end
return