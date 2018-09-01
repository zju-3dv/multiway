function [ L12, E12, nP1, nP2 ] = updateL12E12( cand_matchlist )

featIdx1 = unique(cand_matchlist(:,1));
featIdx2 = unique(cand_matchlist(:,2));
nP1 = length(featIdx1);
nP2 = length(featIdx2);

L12 = zeros(size(cand_matchlist));
for i = 1:length(featIdx1)
    L12(cand_matchlist(:,1) == featIdx1(i), 1) = i;
end
for i = 1:length(featIdx2)
    L12(cand_matchlist(:,2) == featIdx2(i), 2) = i;
end

E12 = zeros(nP1, nP2); 
for i = 1:length(L12), E12(L12(i,1), L12(i,2)) = 1; end;
