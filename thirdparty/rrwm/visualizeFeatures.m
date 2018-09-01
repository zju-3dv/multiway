function visualizeFeatures( frame, varargin )

offset = [0 0];
style = 'point';
colorCode = 'y';
weight = 1;

norMarkerSize = 4;
minMarkerSize = 2; 
maxMarkerSize = 25;
norLineWidth = 2;
minLineWidth = 1;
maxLineWidth = 10;

for k=1:2:length(varargin)
  opt=lower(varargin{k}) ;
  arg=varargin{k+1} ;
  switch opt
    case 'offset'
      offset = arg;
      if numel(offset) == 1
          offset(2) = offset(1);
      end
    case 'markersize'
        norMarkerSize = arg;
    case 'linewidth'
        norLineWidth = arg;    
    case 'style'
      style = arg ;
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
    frame(1,:) = frame(1,:) + offset(1);
    frame(2,:) = frame(2,:) + offset(2);
end

%visScale = 1.0;
%frame(3:6,:) = frame(3:6,:) * visScale;
            
if size(colorCode,2) == 1 && size(weight,2) == 1
    switch style
        case 'point'
            plot(frame(1,:), frame(2,:),'o','MarkerEdgeColor','k',...
                'MarkerFaceColor',colorCode ,'MarkerSize', norMarkerSize);
            
        case 'frame'
            vl_plotframe( frame, 'Color', 'k', 'LineWidth', 4);
            vl_plotframe( frame, 'Color', colorCode, 'LineWidth', 2);

        otherwise
    end
    
else % show each feature with different colors or sizes
    
    minW = min(weight); maxW = max(weight);
    intv_m = (maxMarkerSize-minMarkerSize)/(maxW - minW);
    if (maxW - minW) < 1e-6 
        markerSize = ones(1,size(frame,2)) * norMarkerSize;
        lineWidth = ones(1,size(frame,2)) * norLineWidth;
    else
        markerSize = minMarkerSize+round(intv_m*(weight-minW));
        %lineWidth = minLineWidth+round(intv_l*(weight-minW));
        lineWidth = ones(1,size(frame,2)) * norLineWidth;
    end
    
    if size(colorCode,2) == 1
        colorCode = repmat(colorCode, [1 size(frame,2)]);
    end
    
    switch style
        case 'point'
            for i=1:size(frame,2)
                plot(frame(1,i),frame(2,i),'o','MarkerEdgeColor','k',...
                    'MarkerFaceColor',colorCode(:,i),'MarkerSize', markerSize(i));
            end

        case 'frame'
            for i=1:size(frame,2)
                vl_plotframe( frame(:,i), 'Color', 'k', 'LineWidth', lineWidth(i)+1);
                vl_plotframe( frame(:,i), 'Color', colorCode(:,i), 'LineWidth', lineWidth(i));
            end

        otherwise
    end
            
                    
end


function outColor = colchar2vec(inColor)

  charValues = 'rgbcmywk'.';  % '
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