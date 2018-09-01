function [ adjmat ] = makeAdjMatOfCompleteGraph( nodes )

adjmat = logical(zeros(max(nodes))); 
adjmat(nodes,nodes) = true;
adjmat(1:size(adjmat)+1:end)= false;