function [ featInfo ] = extract_localfeatures_mcho( filePathName, fparam, varargin )
%   extract local features, frames, and their patches
%   modified to use the frame structure 
%   Jan 2013 by Minsu Cho, INRIA - WILLOW / ENS 
%
%   input:   filename and parameters
%   output:  
%   featInfo.frame( :, featIdx ) - x y a b c orientation scale
%   featInfo.type( featIdx )
%   featInfo.desc( :, featIdx)
%   featInfo.patch( :, featIdx)

bVerbose = false ;
minAreaRatio = 0.01^2;
maxAreaRatio = 0.1^2;
bbox_in = [];
thresholdlevel = 1; % 0: low, 1: normal, 2: high

for k=1:2:length(varargin)
  opt=lower(varargin{k}) ;
  arg=varargin{k+1} ;
  switch opt
    case 'verbose'
      bVerbose = arg ;
    case 'minarearatio'
      minAreaRatio = arg ;
    case 'maxarearatio'
      maxAreaRatio = arg ;
    case 'bbox'
      bbox_in = arg ;
    case 'thresholdlevel'
      thresholdlevel = arg;
    otherwise
      error(sprintf('Unknown option ''%s''', opt)) ;
  end
end

featExt = fparam.featExt;
%% Load features from each image
featInfo.img = imread(filePathName);
if size(featInfo.img,3) == 3
    featInfo.img_gray = rgb2gray(featInfo.img);
elseif size(featInfo.img,3) == 1
    featInfo.img_gray = featInfo.img;
    featInfo.img = repmat( featInfo.img, [ 1 1 3 ]);
else
    err([ 'wrong image file!: ' filePathName ]);
end
featInfo.frame = []; featInfo.type = []; featInfo.desc = []; featInfo.patch = [];
%feat = []; typeFeat = []; desc = []; patch = [];

%%%%%%%%%%%%%%%%% Modified code %%%%%%%%%%%%%%%%%%%%%%
if fparam.patchSize < 1
    fparam.patchSize = fparam.patchSize * min(size(featInfo.img_gray));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(bbox_in)
    areaimg = prod(size(featInfo.img_gray));
else
    areaimg = ( bbox_in(3) - bbox_in(1) ) * ( bbox_in(4) - bbox_in(2) );
end

idxFeatExt = find(fparam.bFeatExtUse); % feature indexes to use

for i=idxFeatExt
    strFeatExt = featExt{i};
    if strcmp(strFeatExt(end-2:end),'_vl') % vlfeat detectors
        [tmpFrame tmpDesc ] = loadfeatures_vlfeat( featInfo.img_gray, strFeatExt, thresholdlevel);
        %disp([ filePathName ' ' featExt{i} ': ' num2str(nTmpFeat)]);
        
    else % other feature detectors
        [tmpFeat,tmpDesc ] = loadfeatures_v3( filePathName, strFeatExt, thresholdlevel );
        if fparam.bEstimateOrientation, nMaxOri = fparam.nMaxOri;
        else nMaxOri = 0; end
        
        [featIdx feat_aug tmpOri affMatrix ] = ...
            compute_norm_trans_image( tmpFeat, featInfo.img, fparam.descScale, fparam.patchSize, nMaxOri );
        % reinitialize feat information
        tmpFrame = [ feat_aug(:,1:2)'; affMatrix(:,1)'; affMatrix(:,4)'; affMatrix(:,2)'; affMatrix(:,5)' ];
        %
    end

    % eleminate too small or large features
    area = pi * abs(tmpFrame(3,:).*tmpFrame(6,:) - tmpFrame(4,:).*tmpFrame(5,:)); % area of ellipse
    area = area / fparam.bFeatScale(i);
    idxValid = ( area > minAreaRatio * areaimg ) & ( area < maxAreaRatio * areaimg );
    tmpFrame = tmpFrame(:,idxValid);
    
    %% descriptor extraction
    descriptor = 'sift';
    if size(tmpFrame,2) > 0 && ( isempty(tmpDesc) || size(tmpFrame,2) ~= size(tmpDesc,2) ) 
        [ tmpFrame tmpDesc ] = extractSIFT(featInfo.img_gray, tmpFrame,...
            'contrastInsensitive', fparam.bContrastInsenstive, 'NBP', fparam.nBP,...
            'descScale', fparam.descScale, 'estimateOrientation', fparam.bEstimateOrientation);
        %[ tmpFrame tmpDesc ] = extractSIFT(featInfo.img_gray, tmpFrame, 'estimateOrientation', true);
    end
    
    featInfo.frame = [ featInfo.frame tmpFrame ];
    featInfo.type = [ featInfo.type i*ones(1,size(tmpFrame,2)) ];
    featInfo.desc = [ featInfo.desc tmpDesc ];
%     if bShowFeat
%         figure('Name',[ featExt{i} ': ' filePathName ],'NumberTitle','off')
%         imshow(featInfo.img_gray);
%         hold on; visualizeFeatures(tmpFrame); hold off;
%     end
end

%hold on; vl_plotframe(frame); hold off;
%vl_plotframe( [ frame(1:2,:); mapFromS(frame(3:5,:)) ] , 'Color', 'g');
% apply the margin area
if fparam.marginRatio > 0
    if isempty(bbox_in)
        bbox_in = [ 1 1 size(featInfo.img,2) size(featInfo.img,1) ];
    end
    bbox_in = bbox_in + [ 1 1 -1 -1 ]*ceil(0.5*fparam.marginRatio*min(size(featInfo.img,1),size(featInfo.img,2)));
end

% eliminate features out of a given bounding box
if ~isempty(bbox_in)
    valid_idx = find( (featInfo.frame(1,:) >= bbox_in(1)) & (featInfo.frame(1,:) <= bbox_in(3))...
        & (featInfo.frame(2,:) >= bbox_in(2)) & (featInfo.frame(2,:) <= bbox_in(4)) );
    featInfo.type = featInfo.type( valid_idx );
    featInfo.frame = featInfo.frame( :, valid_idx);
    featInfo.desc = featInfo.desc( :, valid_idx);
    %viewInfo.patch = viewInfo.patch( :, valid_idx);
end

featInfo.fileName = filePathName; 
featInfo.patch = cell(0); 
featInfo.desc_ref = [];
featInfo.nfeature = size(featInfo.desc,2);