function [ featBin rho_bins theta_bins ] = makeFeatBin( attrE, nRho_win, nTheta_win )
% make a polar feature representaion using bins
% output type: single

if nargin < 2
    nRho_win = 3;
    nTheta_win = 3;
end

%if nargin < 4
    %rho_bins = sqrt(3).^([-2:5 Inf]);
    rho_bins = sqrt(2).^([-3:8 Inf]);
    theta_bins = [ -pi + pi/18 : pi/9 : pi - pi/18 ]; % notice the last overlap!
    %theta_bins = [ -pi + pi/36 : pi/18 : pi - pi/36 ]; % notice the last overlap!
%end

%if nargin < 6
    rho_win = gausswin(nRho_win);   % odd number of elements
    theta_win = gausswin(nTheta_win); % odd number of elements
%end

featBin = computePolarBins( attrE, rho_bins', theta_bins', rho_win, theta_win ); 

% nRho_bins = numel(rho_bins);
% nTheta_bins = numel(theta_bins);
% 
% l_tb = floor(numel(theta_win)/2);
% l_rb = floor(numel(rho_win)/2);
% 
% featBin = single(zeros(nTheta_bins+nRho_bins, size(attrE,2)));
% 
% for i=1:size(attrE,2)
%         
%     rho = attrE(1,i);
%     theta = attrE(2,i);
%     %scale = attrE(3,i);
%     if rho == 0 && theta == 0, continue; end
%     % coding
%     i_rb = find(rho < rho_bins, 1, 'first');
%     i_tb = find(theta < [ theta_bins Inf ], 1, 'first');
%     
% 
%     % truncated binning for distance
%     irb_win = i_rb-l_rb:i_rb+l_rb;
%     ivalid = irb_win > 0 & irb_win <= nRho_bins;
%     featBin(irb_win(ivalid),i) = rho_win(ivalid);
% 
%     % circular binning for angle
%     itb_win = i_tb-l_tb:i_tb+l_tb;
%     itb_win(itb_win<1) = itb_win(itb_win<1) + nTheta_bins;
%     itb_win(itb_win>nTheta_bins) = itb_win(itb_win>nTheta_bins) - nTheta_bins;
%     featBin(nRho_bins + itb_win, i) = theta_win(:); 
%         
% end




