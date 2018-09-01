function [ idx ID ] = make_groups( group )
%% Input
% group : ( # match x # groups ) : bianry matrix
%% Output
% idx : ( # match x # groups ) : concatenation of indices of groups elements
% ID : ( # match x # groups ) : concatenation of group ID

nGroup = size(group,2);
idx = [];
ID = [];
for i = 1:nGroup
    tmp = find(group(:,i));
    idx = [idx; tmp];
    ID = [ID; i*ones(size(tmp))];
end