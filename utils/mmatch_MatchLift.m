function [X,Y,info] = mmatch_MatchLift(W,dimGroup,k)

OPT = linear_cons(k,dimGroup);
OPT.stateDims = dimGroup;
OPT.W = sparse(size(W,1)+1,size(W,2)+1);
OPT.W(2:end,2:end) = W;

[X,info] = partial_map_admm_2(OPT);
X = X(2:end,2:end)>0.5;
Y = [];

function [OPT] = linear_cons(m, stateDims)

dim = sum(stateDims)+1;
offsets = ones(1, length(stateDims)+1);

for i = 1:length(stateDims)
    offsets(i+1) = offsets(i) + stateDims(i);
end

% diagonal blocks
ids = 2:dim;
rowsA = 1:(2*dim-1);
colsA = [1, ((ids-1)*dim + ids), ((ids-1)*dim +1)];
valsA = ones(1, 2*dim-1);
b = [m, ones(1, 2*dim-2)];
numEqCons = 2*dim-1;

for i = 1:length(stateDims)
    for j = (offsets(i)+1):offsets(i+1)
        for k = (j+1):offsets(i+1)
            numEqCons = numEqCons + 1;
            id = (k-1)*dim + j;
            colsA = [colsA, id];
            rowsA = [rowsA, numEqCons];
            valsA = [valsA, 1];
            b = [b, 0];
        end
    end

end
OPT.A = sparse(rowsA, colsA, valsA, length(b), dim*dim);
OPT.b = b';