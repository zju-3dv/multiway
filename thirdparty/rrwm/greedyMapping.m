function flag = greedyMapping(match,score,nMax)

% greedy matching 

nMatch = size(match,2);

% identify groups
list = match(1,:);
list_uni = unique(list);
nGroup = length(list_uni);
group = logical(sparse(nGroup,nMatch));
for i=1:nGroup
    group(i,list==list_uni(i)) = true;
end
group1 = group';

list = match(2,:);
list_uni = unique(list);
nGroup = length(list_uni);
group = logical(sparse(nGroup,nMatch));
for i=1:nGroup
    group(i,list==list_uni(i)) = true;
end
group2 = group';

%  greedy selection
flag = false(length(score),1);
count = 0;
[max_value,max_ind] = max(score);
while max_value > -inf && count < nMax
    count = count + 1;
    flag(max_ind) = true;
    score((group1(:,group1(max_ind,:)))) = -inf;
    score((group2(:,group2(max_ind,:)))) = -inf;
    [max_value,max_ind] = max(score);
end
