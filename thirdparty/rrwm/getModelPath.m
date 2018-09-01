function [ path ] = getModelPath( conf, iExp )
if nargin < 2
    iExp = 1;
end
path = sprintf('%s/%s_%s-model%03d_%s_%s.mat',...
    conf.resDir,conf.prefix,conf.class, iExp, conf.referenceType,conf.learningType);

end

