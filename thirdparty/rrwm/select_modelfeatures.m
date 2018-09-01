function [ view_m ] = select_modelfeatures(imgset, nMaxFeat)
% function to select reliable features using GM and NMS
% by Minsu Cho

rNMS = 0.10;          % non-maximal supression radius ratio
rScale = 10.0;          % non-maximal supression scale ratio

set_param_GM;
set_alg;
% extract features
viewInfo1 = extract_localfeatures_mcho( imgset(1).imgpath, fparam, 'bbox', imgset(1).bbox,...
    'minAreaRatio', 0.02, 'maxAreaRatio', 0.2, 'thresholdlevel', 2 );
cdata.view(1) = viewInfo1;  
cdata.fparam = fparam;  cdata.aparam = aparam;  cdata.mparam = mparam;
bbox_in_m = imgset(1).bbox;
if isempty(bbox_in_m), lNMS = rNMS*max(size(viewInfo1.img));
else lNMS = rNMS*max(bbox_in_m(3:4)); end

offset_x = size(cdata.view(1).img, 2);
sV = zeros(1,size(viewInfo1.frame,2)); % score of each node
nV = zeros(1,size(viewInfo1.frame,2)); % hit number
hfig1 = figure('Name','model initialization', 'NumberTitle', 'off'); 
for ic = 2:numel(imgset)
    % extract features
    %viewInfo2 = extract_localfeatures_mcho( imgset(ic).imgpath, fparam, 'bbox', imgset(ic).bbox, 'minAreaRatio', 0.02, 'thresholdlevel', 1  );
    viewInfo2 = extract_localfeatures_mcho( imgset(ic).imgpath, fparam, 'thresholdlevel', 1  );
    cdata.view(2) = viewInfo2;
    bbox_in_t = imgset(ic).bbox;
    
    % visualize it
    imgInput = appendimages( cdata.view(1).img, cdata.view(2).img );
    clf; imshow(rgb2gray(double(imgInput)./255)); hold on;
    visualizeFeatures(viewInfo1.frame, 'style', 'point');
    visualizeFeatures(viewInfo2.frame, 'style', 'point', 'offset', [offset_x 0]);
    box_handle1 = rectangle('position', bbox2rect(bbox_in_m) );
    set(box_handle1, 'edgecolor','y', 'linewidth',5);
    bbox_in_t = bbox_in_t + [ offset_x 0 offset_x 0 ];
    box_handle2 = rectangle('position', bbox2rect(bbox_in_t) );
    set(box_handle2, 'edgecolor','y', 'linewidth',5);

    %% perform matching
    res = wrapper_matchingBasic( alg(1), cdata );
    fprintf('%10s - score:%.5f (%.5f), %.3fsec \n', alg(1).strName, res.score_raw, res.score_GM, res.time);      
    
    %scoreMatch = res.X_raw(find(res.X))'; scoreMatch = scoreMatch./sum(scoreMatch); 
    scoreMatch = res.vSimVal(find(res.X));
    matchSol = res.cand_matchlist(:,find(res.X));
    visualizeMatches( viewInfo1.frame, viewInfo2.frame, matchSol, 'offset', [ offset_x 0 ], 'style', 'point', 'colorcode', 'y', 'weight', scoreMatch );  
    %visualizeFeatures(viewInfo1.frame(:,matchSol(1,:)), 'style', 'point', 'colorcode', 'r');    
    
    % accumulate the GM confidences for each node
    sV(matchSol(1,:)) = sV(matchSol(1,:)) + scoreMatch;
    nV(matchSol(1,:)) = nV(matchSol(1,:)) + 1;
    %visualizeFeatures(viewInfo1.frame, 'style', 'point', 'colorcode', 'm', 'weight', sV);
    sV2 = sV./nV; % consider avg scores
    sV2(nV<0.5*mean(nV)) = 0; % but, exclude nodes with too few hits
    %% select consistent features from the model image
    [ sV_sorted iV_sorted ] = sort(sV2,'descend');
    sel_iV = []; sel_sV = []; nSelFeat = 0; bValid = 0; 
    % compute the radius for non-max suppresion
    for i=1:numel(sV)
        % eliminate redundant matches based on the position
        % accumulate features accoring to the rank of sqdist (ascending)    
        if nSelFeat >= nMaxFeat, break; end
        if i == 1 % the best one
            bValid = 1;
        else    
            frame_cur = viewInfo1.frame(:,iV_sorted(i));    
            frames_sel = viewInfo1.frame(:,sel_iV);
            %finds the points within the search radius
            cand_iidx = BruteSearchMex(frames_sel(1:2,:),frame_cur(1:2),'r',lNMS);
            if isempty(cand_iidx) 
                bValid = 1; % if there's no feature with close distance
            else
                scales_change = getFrameScale(frames_sel(:,cand_iidx))/getFrameScale(frame_cur);
                if all( (scales_change > rScale) | (scales_change < 1/rScale) )
                    bValid = 1; % if there's no features with close scale
                    %disp('scale!');
                end
            end
        end

        if bValid
            nSelFeat = nSelFeat + 1;
            % insert it into the selected xy list
            sel_iV(nSelFeat) = iV_sorted(i);
            sel_sV(nSelFeat) = sV_sorted(i);
            visualizeFeatures( viewInfo1.frame(:,sel_iV(nSelFeat)), 'style', 'frame', 'colorcode', 'r' );
            %visualizeFeatures( viewInfo1.frame(:,sel_iV(nSelFeat)), 'style', 'point', 'colorcode', 'r' );
            text( viewInfo1.frame(1,sel_iV(nSelFeat)), viewInfo1.frame(2,sel_iV(nSelFeat)),...
                num2str(nSelFeat), 'FontSize',20); 
            %pause;
            bValid = 0;
        end
    end
    pause(0.5);
    %sel_sV
    %pause;
end

% show the image
%figure('Name','Select nodes...', 'NumberTitle', 'off'); 
clf;    imshow(rgb2gray(viewInfo1.img)); hold on;
visualizeFeatures( viewInfo1.frame, 'style', 'point');
box_handle1 = rectangle('position', bbox2rect(bbox_in_m));
set(box_handle1, 'edgecolor','y', 'linewidth',5);
visualizeFeatures( viewInfo1.frame(:,sel_iV), 'style', 'frame', 'colorcode', 'r' );
visualizeFeatures( viewInfo1.frame(:,sel_iV), 'style', 'point', 'colorcode', 'r' );
for i=1:numel(sel_iV)
    text( viewInfo1.frame(1,sel_iV(i)), viewInfo1.frame(2,sel_iV(i)), num2str(i), 'FontSize',20); 
end

pause(0.5);
close(hfig1);

view_m = viewInfo1;
view_m.type = view_m.type( sel_iV );
view_m.frame = view_m.frame(:, sel_iV);
view_m.desc = view_m.desc(:, sel_iV);

%% select a bouding box
%rect = getrect;
rect = bbox2rect(bbox_in_m);
nGrid = 1;
x_min = rect(1);    y_min = rect(2);
x_intv = round(rect(3)/nGrid);  y_intv = round(rect(4)/nGrid);
x_max = x_min + x_intv*nGrid;   y_max = y_min + y_intv*nGrid;

% Display the subsetted image with appropriate axis ratio
%figure; 
%image(view_m.img(rect(2):rect(2)+rect(4), rect(1):rect(1)+rect(3),:)); axis image
%pause;
bValid_x = view_m.frame(1,:) >= rect(1) & view_m.frame(1,:) < rect(1)+rect(3);
bValid_y = view_m.frame(2,:) >= rect(2) & view_m.frame(2,:) < rect(2)+rect(4);
bValid = bValid_x & bValid_y; 

%plot(view_m.feat(find(bValid),1),view_m.feat(find(bValid),2),'o','MarkerEdgeColor','k',...
%    'MarkerFaceColor','r','MarkerSize',10);
idx_selfeat = find(bValid);

view_m.type = view_m.type(idx_selfeat);
view_m.frame = view_m.frame(:, idx_selfeat);
view_m.desc = view_m.desc(:, idx_selfeat);