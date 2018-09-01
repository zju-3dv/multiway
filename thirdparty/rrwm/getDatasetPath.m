function [ path ] = getDatasetPath( conf )

path = sprintf('%s/%s_%s_dataset.mat',...
    conf.resDir,conf.prefix, conf.class);

end
