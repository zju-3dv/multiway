function [ edges ] = makeStarGraphEdges( nodes, frame )

center_pt = mean(frame(1:2, nodes),2);
d = sum((frame(1:2,:) - repmat(center_pt,1,numel(nodes))).^2,1);
[ tmp iv ] = min(d);
edges = [iv*ones(1,numel(nodes)-1); 1:iv-1, iv+1:numel(nodes)]; 