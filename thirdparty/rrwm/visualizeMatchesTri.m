function visualizeMatchesTri( frame1, frame2, match, varargin )


offset = [0 0];
style = 'point';
colorCode = 'y';
weight = 1;

norMarkerSize = 12;
minMarkerSize = 2; 
maxMarkerSize = 30;
norLineWidth = 2;
minLineWidth = 1;
maxLineWidth = 10;
indicator = ones(1,size(match,2));

for k=1:2:length(varargin)
  opt=lower(varargin{k}) ;
  arg=varargin{k+1} ;
  switch opt
    case 'offset'
      offset = arg;
      if numel(offset) == 1
          offset(2) = offset(1);
      end
    case 'style'
      style = arg ;
    case 'indicator'
      indicator = arg ;  
    case 'colorcode'
      colorCode = arg ;  
    case 'weight'
      weight = arg ;    
    otherwise
      error(sprintf('Unknown option ''%s''', opt)) ;
  end
end

if ischar(colorCode)
    colorCode = colchar2vec(colorCode)';
end

if nnz(offset) > 0
    frame2(1,:) = frame2(1,:) + offset(1);
    frame2(2,:) = frame2(2,:) + offset(2);
end


%visScale = 6.0;

minW = min(weight); maxW = max(weight);
intv_m = (maxMarkerSize-minMarkerSize)/(maxW - minW+eps);
intv_l = (maxLineWidth-minLineWidth)/(maxW - minW+eps);

if minW == maxW
    markerSize = ones(1,size(match,2)) * norMarkerSize;
    lineWidth = ones(1,size(match,2)) * norLineWidth;
else
    markerSize = minMarkerSize+round(intv_m*(weight-minW));
    %lineWidth = minLineWidth+round(intv_l*(weight-minW));
    lineWidth = ones(1,size(match,2)) * norLineWidth;
end

tri = delaunay(frame1(1,:),frame1(2,:));
%triplot(tri,frame1(1,:),frame1(2,:),'color','w','linewidth',4);
triplot(tri,frame1(1,:),frame1(2,:),'color','r','linewidth',4);

%triplot(tri,frame2(1,match(1,:)),frame2(2,match(2,:)),'color','w','linewidth',4);
triplot(tri,frame2(1,match(2,:)),frame2(2,match(2,:)),'color','r','linewidth',4);

idxFalse = find(indicator == 0);
for j=1:numel(idxFalse)
    k = idxFalse(j);
    % find connected nodes to k
    tmp = find( tri(:,1) == k | tri(:,2) == k | tri(:,3) == k );
    tmp2 = tri(tmp,:);
    idxConnec = setdiff( unique(tmp2(:)), k );

    for i=1:numel(idxConnec)
        line( [ frame1(1,k) frame1(1,idxConnec(i))],...
            [frame1(2,k) frame1(2,idxConnec(i))], 'color','w','linewidth',4);
        line( [ frame1(1,k) frame1(1,idxConnec(i))],...
            [frame1(2,k) frame1(2,idxConnec(i))],'linestyle',':', 'color','k','linewidth',4);
        line( [ frame2(1,match(2,k)) frame2(1,match(2,idxConnec(i)))] ,...
            [frame2(2,match(2,k)) frame2(2,match(2,idxConnec(i)))], 'color','w','linewidth',4);
        line( [ frame2(1,match(2,k)) frame2(1,match(2,idxConnec(i)))] ,...
            [frame2(2,match(2,k)) frame2(2,match(2,idxConnec(i)))],'linestyle',':', 'color','k','linewidth',4);
    end
end


for i = 1:size(match,2)
    if size(colorCode,2) < i
        col = rand(3,1);
    else
        col = colorCode(:,i);
    end
    
    switch style
        case 'point'
            plot([ frame1(1,match(1,i)) frame2(1,match(2,i))],...
                [ frame1(2,match(1,i)) frame2(2,match(2,i))],...
                'o','Color', col, 'MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize', markerSize(i));
            
        case 'frame'
            vl_plotframe( frame1(:,match(1,i)), 'Color', 'k', 'LineWidth', lineWidth(i)+1);
            vl_plotframe( frame1(:,match(1,i)), 'Color', col, 'LineWidth', lineWidth(i));
            vl_plotframe( frame2(:,match(2,i)), 'Color', 'k', 'LineWidth', lineWidth(i)+1);
            vl_plotframe( frame2(:,match(2,i)), 'Color', col, 'LineWidth', lineWidth(i));
            %plot([ frame1(1,match(1,i)) frame2(1,match(2,i))],...
            %    [ frame1(2,match(1,i)) frame2(2,match(2,i))],...
            %    '-o','Color', col, 'MarkerEdgeColor','k','MarkerFaceColor',col,'MarkerSize', markerSize(i));
    
        otherwise
    end
            
end        

for j=1:numel(idxFalse)
    k = idxFalse(j);
    plot([ frame1(1,match(1,k)) frame2(1,match(2,k))],...
                [ frame1(2,match(1,k)) frame2(2,match(2,k))],...
                'o','Color', 'k', 'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize', markerSize(k));
end

function outColor = colchar2vec(inColor)

  charValues = 'rgbcmywk'.';  % #'
  rgbValues = [eye(3); 1-eye(3); 1 1 1; 0 0 0];
  assert(~isempty(inColor),'convert_color:badInputSize',...
         'Input argument must not be empty.');

  if ischar(inColor)  % input is a character string

    [isColor,colorIndex] = ismember(inColor(:),charValues);
    assert(all(isColor),'convert_color:badInputContents',...
           'String input can only contain the characters ''rgbcmywk''.');
    outColor = rgbValues(colorIndex,:);
  
  else  % input is an invalid type

    error('convert_color:badInputType',...
          'Input must be a character or numeric array.');

  end
