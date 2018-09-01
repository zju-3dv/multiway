function [ edges ] = makeDelaunayGraphEdges( nodes, frame )
% working on....
edges = [];

delaunayTRI = DelaunayTri(frame(1,:),frame(2,:));
%triplot(delaunayTRI,'r-','LineWidth',3);
edges = edges(delaunayTRI);
edges = edges';