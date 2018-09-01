function drawWeightedGraph( frame_m, offset, color_m, edges_m, wV, wE )
hold on;
% show node shapes (ellipses)
%visualizeFeatures( frame_m, 'offset', offset, 'style', 'frame', 'colorcode', color_m);

% visualize edges with their weights
visualizeEdges( frame_m, edges_m, 'offset', offset, 'weight', wE, 'linestyle', '-', 'colorcode', 'b');

% visualize nodes with their weights
visualizeFeatures( frame_m, 'offset', offset, 'style', 'point', 'colorcode', color_m, 'weight', wV, 'markersize', 15);

for i = 1:size(frame_m,2)    
    %text(frame_m(1,i) + offset(1) , frame_m(2,i) + offset(2) , num2str(i), 'color','k','HorizontalAlignment', 'center','FontSize',15);        
end

hold off;