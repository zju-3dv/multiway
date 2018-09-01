function flag = greedyMatch(match,score,nMax)

if nargin < 3
    nMax = inf;
end

%  greedy selection
flag = false(length(score),1);
count = 0;
[max_value,max_ind] = max(score);
while max_value > -inf && count < nMax
    count = count + 1;
    flag(max_ind) = true;
    score(match(1,:)==match(1,max_ind)) = -inf;
    score(match(2,:)==match(2,max_ind)) = -inf;
    [max_value,max_ind] = max(score);
end
