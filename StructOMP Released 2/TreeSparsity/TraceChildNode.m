function Index=TraceChildNode(Ind, Edge4);

if size(Ind, 2)~=1
    Ind=Ind';
end
[m, d]=size(Edge4);

num=length(Ind);

stop=1;
Index=unique(Ind);
while stop
    Index0=Index;
    Temp=Edge4(Index, 3:end);
    Index=[Index; Temp(:)];
    Index(find(Index==0))=[];
    Index=unique(Index);

    if length(Index0)==length(Index) 
        if Index0==Index
            break
        end
    end
end