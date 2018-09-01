function [l,p] = IsInBox(p,bbox)

l = p(1,:) >= bbox(1) & p(1,:) <= bbox(3) & ...
    p(2,:) >= bbox(2) & p(2,:) <= bbox(4);

p(1,:) = min(max(p(1,:),bbox(1)),bbox(3));
p(2,:) = min(max(p(2,:),bbox(2)),bbox(4));