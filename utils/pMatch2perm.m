function M = pMatch2perm(pMatch)

nFeature = zeros(size(pMatch,1),1);
filename = cell(size(pMatch,1),1);
for i = 1:size(pMatch,1)
    for j = i+1:size(pMatch,2)
        if ~isempty(pMatch(i,j).nFeature)
            nFeature(i) = pMatch(i,j).nFeature(1);
            nFeature(j) = pMatch(i,j).nFeature(2);
        end
        if ~isempty(pMatch(i,j).filename)
            filename(i) = pMatch(i,j).filename(1);
            filename(j) = pMatch(i,j).filename(2);
        end
    end
end
cumIndex = cumsum([0; nFeature]);

M = sparse(cumIndex(end),cumIndex(end));
for i = 1:size(pMatch,1)
    for j = i+1:size(pMatch,2)
        if ~isempty(pMatch(i,j).matchInfo)
            matchList = double(pMatch(i,j).matchInfo.match);
            M(cumIndex(i)+1:cumIndex(i+1),cumIndex(j)+1:cumIndex(j+1)) = ...
                sparse(matchList(1,:),matchList(2,:),pMatch(i,j).X,nFeature(i),nFeature(j));
        end
    end
end
M = M + M';