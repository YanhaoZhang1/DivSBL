function west=OverlapLasso(X, y, varGroupMatrix, lambda);

[nVars, nGroups]=size(varGroupMatrix);

vInd = find(varGroupMatrix==1);
nSubVars = length(vInd);
v = zeros(nSubVars,1);

alpha = zeros(nGroups,1);
funObj = @(w)SquaredError(w,X,y);
penalizedFunObj = @(vAlpha)overLassoLoss(vAlpha,varGroupMatrix,lambda,funObj);

% Set up sub-variable groups
subGroups = zeros(nSubVars,1);
offset = 0;
for g = 1:nGroups
    subGroupLength = sum(varGroupMatrix(:,g));
    subGroups(offset+1:offset+subGroupLength) = g;
    offset = offset+subGroupLength;
end

% Set up L_1,inf Projection Function
[groupStart,groupPtr] = groupl1_makeGroupPointers(subGroups);
funProj = @(vAlpha)auxGroupL2Project(vAlpha,nSubVars,groupStart,groupPtr);

% Solve with PQN
% fprintf('\nComputing over-lasso regularized linear regression parameters...\n');
options.verbose=0; options.optTol=1e-6;
vAlpha = minConF_PQN(penalizedFunObj,[v;alpha],funProj, options);

% Extract parameters from augmented vector
v = vAlpha(1:nSubVars);
v(abs(v) < 1e-4) = 0;

% Form sub-weight matrix vFull, and weight vector w
vFull = zeros(nVars,nGroups);
vFull(vInd) = v;
west = sum(vFull,2);
return
