function [frame,desc] = loadfeatures_vlfeat(img_gray, detector, thresholdlevel)

% Load local features using vlfeat library
% Author: Minsu Cho
% Modified by Xiaowei Zhou

param.EstimateAffineShape = false;
param.EstimateOrientation = false;

switch lower(detector)
    
    case 'dog_vl' % DOG features
        if thresholdlevel == 2
            param.DoG_PeakThreshold = 1.5;
        else
            param.DoG_PeakThreshold = 0.5;
        end
        frame = vl_covdet(single(img_gray), 'Method', 'DoG',...
            'EstimateAffineShape', param.EstimateAffineShape,...
            'EstimateOrientation', param.EstimateOrientation,...
            'PeakThreshold', param.DoG_PeakThreshold);%,...'verbose') ;
          
    case 'hes_vl'
        if thresholdlevel == 2
            param.Hessian_PeakThreshold = 25.0; %?
        else
            param.Hessian_PeakThreshold = 2.0;
        end
        frame = vl_covdet(single(img_gray), 'Method', 'Hessian',...
            'EstimateAffineShape', param.EstimateAffineShape,...
            'EstimateOrientation', param.EstimateOrientation,...
            'PeakThreshold', param.Hessian_PeakThreshold);%, 'verbose');
        
    case 'heslap_vl'
        if thresholdlevel == 2
            param.HessianLaplace_PeakThreshold = 130.0; %?
        else
            param.HessianLaplace_PeakThreshold = 30.0;
        end
        frame = vl_covdet(single(img_gray), 'Method', 'HessianLaplace',...
            'EstimateAffineShape', param.EstimateAffineShape,...
            'EstimateOrientation', param.EstimateOrientation,...
            'PeakThreshold', param.HessianLaplace_PeakThreshold, 'verbose');
        
    case 'harlap_vl'
        if thresholdlevel == 2
            param.HarrisLaplace_PeakThreshold = 5.0; %?
        else
            param.HarrisLaplace_PeakThreshold = 0.0;
        end
        frame = vl_covdet(single(img_gray), 'Method', 'HarrisLaplace',...
            'EstimateAffineShape', param.EstimateAffineShape,...
            'EstimateOrientation', param.EstimateOrientation,...
            'PeakThreshold', param.HarrisLaplace_PeakThreshold );
        
    case 'mser_vl'
        [~,tframes] = vl_mser(uint8(img_gray)) ;
        tframes = vl_ertr(tframes);
        K = size(tframes,2);
        % transform the format to the format of translation 2 +affine 4
        frame = zeros(6,K) ;
        frame(1:2,:) = tframes(1:2,:);
        frame(3:6,:) = mapFromS(tframes(3:5,:));
   
    otherwise     
        error('Unknown feature type "%s"',detector);
end

desc = [];
        
end


% --------------------------------------------------------------------
function A = mapFromS(S)
% --------------------------------------------------------------------
% Returns the (stacking of the) 2x2 matrix A that maps the unit circle
% into the ellipses satisfying the equation x' inv(S) x = 1. Here S
% is a stacked covariance matrix, with elements S11, S12 and S22.

tmp = sqrt(S(3,:)) + eps ;
A(1,:) = sqrt(S(1,:).*S(3,:) - S(2,:).^2) ./ tmp ;
A(2,:) = zeros(1,length(tmp));
A(3,:) = S(2,:) ./ tmp ;
A(4,:) = tmp ;

end
