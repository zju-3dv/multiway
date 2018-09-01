function [group1,group2] = make_group12(match)
% make group1 and group2 based on the match list
% group1 = ( # of matches ) by ( # of unique features in view 1 )
% group2 = ( # of matches ) by ( # of unique features in view 2 )

nMatch = size(match,2);

featList = match(1,:);
featList_unique = unique(featList);
nGroup = length(featList_unique);
group = logical(sparse(nGroup,nMatch));
for i=1:nGroup
    group(i, featList == featList_unique(i)) = true;
end
group1 = group';

featList = match(2,:);
featList_unique = unique(featList);
nGroup = length(featList_unique);
group = logical(sparse(nGroup,nMatch));
for i=1:nGroup
    group(i, featList == featList_unique(i)) = true;
end
group2 = group';