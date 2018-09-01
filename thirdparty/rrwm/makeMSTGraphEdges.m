function [ edges ] = makeMSTGraphEdges( nodes, frame )
% working on....
edges = [];

% edges by minimum spanning tree
% PV = [];
% for i=1:numel(nodes)-1
%     for j=i+1:numel(nodes)
%         d = sqrt(sum((frame(1:2,i) - frame(1:2,j)).^2));
%         PV = [ PV; i j d ];
%     end
% end
% [w T] = kruskal(PV);
% [ edges(:,1) edges(:,2) ] = find(triu(T,1));

edges = edges';