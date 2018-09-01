function [ featInfo ] = extract_localfeatures_pts( x , varargin )
%   extract local features, frames, and their patches
%   modified to use the frame structure 
%   Jan 2013 by Minsu Cho, INRIA - WILLOW / ENS 
%
%   input:   x
%   output:  
%   featInfo.frame( :, featIdx ) - x y a b c orientation scale
%   featInfo.type( featIdx )
%   featInfo.desc( :, featIdx)

bVerbose = false ;
bbox_in = [];
scale = 0.1;

for k=1:2:length(varargin)
  opt=lower(varargin{k}) ;
  arg=varargin{k+1} ;
  switch opt
    case 'verbose'
      bVerbose = arg ;
    case 'bbox'
      bbox_in = arg ;
    case 'scale'
      scale = arg;
    otherwise
      error(sprintf('Unknown option ''%s''', opt)) ;
  end
end

featInfo.frame = []; featInfo.type = []; featInfo.desc = []; featInfo.patch = [];

tmpFrame = frame2oell([ x; scale*ones(1,size(x,2))]);

featInfo.frame = tmpFrame;
featInfo.type = ones(1,size(tmpFrame,2));
featInfo.desc = zeros(1,size(tmpFrame,2));

if 0
    figure('Name', 'Pts','NumberTitle','off');
    hold on; visualizeFeatures(tmpFrame,'style','frame'); hold off;
    pause;
end

