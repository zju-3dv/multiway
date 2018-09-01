function [ path ] = getResultPath( conf )

path = sprintf('%s/%s_%s-res_%s_%s.mat',...
    conf.resDir,conf.prefix, conf.class,  conf.referenceType,conf.learningType);

end
