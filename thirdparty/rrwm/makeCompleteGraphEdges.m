function [ edges ] = makeCompleteGraphEdges( nodes )

edges = [];
for i=1:numel(nodes)-1
    for j=i+1:numel(nodes)
        edges = [ edges; i j ];
    end
end
edges = edges';