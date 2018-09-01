function showMatchingResult( cdata, res, colorCode )
% visualize it
if nargin < 3
    colorCode = makeColorCode(nnz(res.X));
end

imshow(rgb2gray(double( appendimages( cdata.view(1).img, cdata.view(2).img ))./255)); hold on;
offset_x = size(cdata.view(1).img,2);
%visualizeFeatures(cdata.view(1).frame, 'style', 'point', 'colorcode', 'g');
%visualizeFeatures(cdata.view(2).frame, 'style', 'point', 'colorcode', 'g', 'offset', [offset_x 0] );

% box_handle1 = rectangle('position', bbox2rect(bbox_in1) );
% set(box_handle1, 'edgecolor','y', 'linewidth',5);
% bbox_in2 = bbox_in2 + [ 1 0 1 0 ]*size(cdata.view(1).img, 2);
% box_handle2 = rectangle('position', bbox2rect(bbox_in2) );
% set(box_handle2, 'edgecolor','y', 'linewidth',5);

visualizeMatchesTri( cdata.view(1).frame, cdata.view(2).frame, res.cand_matchlist(:,find(res.X)), ...
    'offset', [size(cdata.view(1).img,2) 0], 'style', 'point',...
    'colorcode', colorCode, 'weight', res.X_raw(find(res.X)), 'indicator', res.indicator ); 

%visualizeMatches( cdata.view(1).frame, cdata.view(2).frame, res.cand_matchlist(:,find(res.X)), ...
%    'offset', [size(cdata.view(1).img,2) 0], 'style', 'frame', 'colorcode', colorCode, 'weight', res.X_raw(find(res.X)) ); 
%visualizeSIFT2(svm.wV, cdata.view(1).frame);

%pause;

%% select consistent features from the model image
% idxMatch = find(res.X);
% scoreMatch = res.X_raw(idxMatch);
% [ tmp iisorted ] = sort(scoreMatch,'descend');
% %iisorted = iisorted(1:10);
% %tmp
% matchSol = res.cand_matchlist(1:2, idxMatch(iisorted));
% 
% idxFeat1 = matchSol(1,:); idxFeat2 = matchSol(2,:);
% frame1 = cdata.view(1).frame(:,idxFeat1);
% frame2 = cdata.view(2).frame(:,idxFeat2);
% frame2(1,:) = frame2(1,:) + size(cdata.view(1).img,2);

% estimate bbox and show it
% extract corresponding features
%matchSol = res.cand_matchlist(:,find(res.X));
%frame2 = cdata.view(2).frame(:,matchSol(2,:));
%close(hFigRes);

end