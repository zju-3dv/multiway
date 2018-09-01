function scale = getFrameScale( frame )
% get scale mutipliers w.r.t unit

scale = abs(frame(3,:).*frame(6,:) - frame(4,:).*frame(5,:));

