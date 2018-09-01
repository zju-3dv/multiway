function [ path ] = getExpPath( conf, postfix )

path = sprintf('%s/%s_%s_exp_%s.mat',...
    conf.resDir,conf.prefix, conf.class, postfix);

end
