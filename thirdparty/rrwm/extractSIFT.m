function [ tframe desc ] = extractSIFT( img_gray, frame, varargin )
% vl_covdet function scales up 6 times to extract descriptors
% To modify it to 2 times, feature scales are changed before extraction.
descScale = 3.0;
estimateOrientation = false;
bContrastInsensitive = true;
NBO = 8;
NBP = 4;

for k=1:2:length(varargin)
  opt=lower(varargin{k}) ;
  arg=varargin{k+1} ;
  switch opt
    case 'descscale'
      descScale = arg;
    case 'estimateorientation'
      estimateOrientation = arg;
    case 'contrastinsensitive'
      bContrastInsensitive = arg;  
    case 'nbo'
      NBO = arg;
    case 'nbp'
      NBP = arg;  
    otherwise
      error(sprintf('Unknown option ''%s''', opt)) ;
  end
end

if size(img_gray,3) > 1
    img_gray = rgb2gray(img_gray);
end

if size(frame,1) ~= 6
    frame = frame2oell(frame);
end

%frame(3:6,:) = frame(3:6,:) * (descScale / 6.0);
%try

if 1
    %[tframe desc] = vl_covdet(single(img_gray), 'frames', frame, 'NBP', NBP, ...
    [tframe desc] = vl_covdet(single(img_gray), 'frames', frame, ...
        'PatchResolution', 50, 'PatchRelativeExtent', descScale, 'PatchRelativeSmoothing', 0.1, ...
      'estimateOrientation', estimateOrientation) ;
%    if bContrastInsensitive
%        HNBO = NBO / 2;
%        for i=1:NBP*NBP
%             desc_ci((i-1)*HNBO+1:i*HNBO,:) = desc((i-1)*NBO+1:(i-1)*NBO+HNBO,:) + desc((i-1)*NBO+HNBO+1:i*NBO,:);
%        end
%        norms = 1./sqrt(sum(desc_ci.^2,1));
%        desc = desc_ci.*repmat(norms,size(desc_ci,1),1);
%    end
else
    [tframe patches] = vl_covdet(single(img_gray), 'descriptor', 'patch', 'frames', frame, ... 
        'PatchResolution', 22, 'PatchRelativeExtent', descScale, 'estimateOrientation', estimateOrientation) ;
    w = sqrt(size(patches,1));
    patches2 = reshape(patches, w,w,[]);
    for i=1:size(patches2,3)
        hog = vl_hog(patches2(:,:,i), 9, 'numOrientations', 8) ;
        desc(:,i) = hog(:)./sqrt(sum(hog(:).^2));
    end
    %w% 22 9, 38 11
    %size(hog)
    %vl_imarraysc(patches2(:,:,1:200)) ;
    %colormap gray;
    %pause;
end

%tframe(3:6,:) = tframe(3:6,:) * (6.0 / descScale);
%frame = [ frame(1:2,:); frame(3:6,:) * 3.0 ];
%featInfo.desc = single(featInfo.desc); % ( this case, already single )                                
if size(tframe,2) ~= size(desc,2) 
    error('feat dim is not consistent with desc dim!');
end