function [BC, BCm]=GetBlocksConnectionMatrix(Bm, Edge4);
[m, maxnum]=size(Bm);

for i=1:m,
    if i/1000==round(i/1000)
        fprintf('i=%d\n', i);
    end
    bi=Bm(i,:); bi(find(bi==0))=[];
    bicon=Edge4(bi,:);    
    bicon=unique(bicon(:));bicon=bicon';
    bicon(find(bicon==0))=[];
    
    tag=ismember(Bm, bicon);
    tag=sum(tag, 2);
    tag(i,1)=0;
    index=find(tag~=0);    
    bci=index;
    BC{i}=bci;   
%     BCm(i,:)=bci;
end

maxnum=0;
for i=1:m,
    tmp=length(BC{i});
    maxnum=max(maxnum, tmp);
end
BCm=zeros(m, maxnum);
for i=1:m,
    tmp=BC{i};
    BCm(i,1:length(tmp))=tmp;
    clear tmp;
end    
return