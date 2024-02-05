function [B, Bmatrix]=CreatTreeBlock(H, W, IdxParent, Step);
%%%% how deep the ancester you want?
m=length(IdxParent);
for i=1:m,
    curr_ind=i;
    list=[i];
    in=1;
    while 1
        curr_ind=IdxParent(list(in));
        if curr_ind<1 | in>Step
            break
        end
        list=[list, curr_ind];
        in=in+1;
    end         
    B{i}=list;
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