function varGroupMatrix=CreatGroupMatrix(H, W, flag);
%%%%% flag==1  1d line structure; each group includes one node and its
%%%%% neighborhood nodes (two)


if flag==1  %%% W=1 in this case
    nGroups=H;nVars = H;
    varGroupMatrix = sparse(nVars,nGroups);
    for g = 1:nGroups
        start=max(1, g-1); stop=min(H, g+1);
        varGroupMatrix(start:stop, g) = ones(stop-start+1, 1);
    end
else
end

