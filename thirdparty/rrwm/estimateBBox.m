function [ bbox ] = estimateBBox(model, frame)
% estimate the bounding box from the matching features
% a simple version

xmin_t = min(frame(1,:));   xmax_t = max(frame(1,:));
ymin_t = min(frame(2,:));   ymax_t = max(frame(2,:));
c_pt = [ xmin_t+xmax_t ymin_t+ymax_t ]*0.5;

% specify the minimal bounding box [xmin ymin width height]
%bbox = [ xmin_t ymin_t xmax_t-xmin_t ymax_t-ymin_t ];
%bbox(1) = bbox(1) + size(cdata.view(1).img,2);
%box_handle = rectangle('position', bbox);
%set(box_handle, 'edgecolor','r', 'linewidth',3);

% inferr bbox using least square fitting
xy_min = model.bbox_affMatrix*[xmin_t-c_pt(1) ymin_t-c_pt(2) 1]';
xy_max = model.bbox_affMatrix*[xmax_t-c_pt(1) ymax_t-c_pt(2) 1]';
xy_min = xy_min + c_pt';   xy_max = xy_max + c_pt';

bbox = [ xy_min(1) xy_min(2) xy_max(1)-xy_min(1) xy_max(2)-xy_min(2) ];
