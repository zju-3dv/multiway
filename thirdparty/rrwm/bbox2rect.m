function [ rect ] = bbox2rect( bbox )
% x1 y1 x2 y2 -> x1 y1 x2-x1 y2-y1 

rect = [bbox(1), bbox(2), bbox(3)-bbox(1), bbox(4)-bbox(2)];
