function visualizeEdges( frame, adjmat, varargin )
% visualze edges connecting feature nodes
% by Minsu Cho

offset = [0 0];
linestyle = '-';
colorCode = 'y';
weight = double(adjmat(:));
vis_weight = 'color';

norLineWidth = 3;
minLineWidth = 1;
maxLineWidth = 10;

if nargin < 2
    error('this function requires two arguments at least!') ;
end

for k=1:2:length(varargin)
  opt=lower(varargin{k}) ;
  arg=varargin{k+1} ;
  switch opt
    case 'offset'
      offset = arg;
      if numel(offset) == 1
          offset(2) = offset(1);
      end
    case 'linestyle'
      linestyle = arg ;
    case 'colorcode'
      colorCode = arg ;  
    case 'weight'
      weight = arg ;
    case 'vis_weight'
      vis_weight = arg;      
    otherwise
      error(sprintf('Unknown option ''%s''', opt)) ;
  end
end

if ischar(colorCode)
    colorCode = colchar2vec(colorCode)';
end

minW = min(abs(weight(adjmat(:)))); maxW = max(abs(weight(adjmat(:))));

if ( maxW - minW ) < 1e-6 || strcmp(vis_weight, 'color')
    lineWidth = double(adjmat(:));
    lineWidth(adjmat) = norLineWidth;
elseif strcmp(vis_weight, 'width')
    %intv_l = (maxLineWidth-minLineWidth)/(maxW - minW);
    intv_l = (maxLineWidth-minLineWidth)/maxW;
    lineWidth = minLineWidth+round(intv_l*abs(weight));
end

if nnz(offset) > 0
    frame(1,:) = frame(1,:) + offset(1);
    frame(2,:) = frame(2,:) + offset(2);
end

% non-positive weight edges
[ i_n, v_n ] = find( adjmat(:) & weight(:) <= 0 );
[ r_n, c_n ] = ind2sub(size(adjmat),i_n);

for i = 1:numel(i_n)
    plot([ frame(1,r_n(i)) frame(1,c_n(i))],...
        [ frame(2,r_n(i)) frame(2,c_n(i))],...
        linestyle,'Color', 'k', 'LineWidth', lineWidth(i_n(i)));
end
% positive weight edges
[ i_p ] = find( adjmat(:) & weight(:) > 0 );
[ r_p, c_p ] = ind2sub(size(adjmat),i_p);

[ wval widx ] = sort(weight(i_p),'ascend');

if strcmp(vis_weight, 'color') && ( maxW - minW ) > 1e-6
    %cmap = colormap('jet');
    cmap = colormap('gray');
    cmap = cmap(end:-1:1,:);
    colormap('gray');
    for i = 1:numel(widx)
        colorCode = cmap( ceil( (wval(i)-minW) * length(cmap) / (maxW-minW) + eps), :);
        plot([ frame(1,r_p(widx(i))) frame(1,c_p(widx(i)))],...
            [ frame(2,r_p(widx(i))) frame(2,c_p(widx(i)))],...
            linestyle,'Color', colorCode, 'LineWidth', lineWidth(i_p(widx(i))));
    end        
else
    for i = 1:numel(widx)
        plot([ frame(1,r_p(widx(i))) frame(1,c_p(widx(i)))],...
            [ frame(2,r_p(widx(i))) frame(2,c_p(widx(i)))],...
            linestyle,'Color', colorCode, 'LineWidth', lineWidth(i_p(widx(i))));
    end        
end

function outColor = colchar2vec(inColor)

  charValues = 'rgbcmywk'.';  %#'
  rgbValues = [eye(3); 1-eye(3); 1 1 1; 0 0 0];
  assert(~isempty(inColor),'convert_color:badInputSize',...
         'Input argument must not be empty.');

  if ischar(inColor)  %# Input is a character string

    [isColor,colorIndex] = ismember(inColor(:),charValues);
    assert(all(isColor),'convert_color:badInputContents',...
           'String input can only contain the characters ''rgbcmywk''.');
    outColor = rgbValues(colorIndex,:);
  
  else  %# Input is an invalid type

    error('convert_color:badInputType',...
          'Input must be a character or numeric array.');

  end
