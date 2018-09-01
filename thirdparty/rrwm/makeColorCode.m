function [ priorColorCode ] = makeColorCode( nCol )

priorColorCode(1,:) = [ 1 0 0 ]; % R
priorColorCode(2,:) = [ 0 1 0 ]; % G
priorColorCode(3,:) = [ 0 0 1 ]; % B
priorColorCode(4,:) = [ 0 1 1 ]; % C
priorColorCode(5,:) = [ 1 0 1 ]; % M
priorColorCode(6,:) = [ 1 1 0 ]; % Y
priorColorCode(7,:) = [ 1 0.5 0 ]; % Y
priorColorCode(8,:) = [ 1 0 0.5 ]; % Y
priorColorCode(9,:) = [ 1 0.5 0.5 ]; % Y
priorColorCode(10,:) = [ 0.5 1 0 ]; % Y
priorColorCode(11,:) = [ 0 1 0.5 ]; % Y
priorColorCode(12,:) = [ 0.5 1 0.5 ]; % Y

nMore = nCol - size(priorColorCode,1);
if nMore > 0 
    priorColorCode(size(priorColorCode,1)+1:nCol,:) = rand(nMore, 3);
end

priorColorCode = priorColorCode';