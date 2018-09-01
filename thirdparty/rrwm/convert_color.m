function outColor = convert_color(inColor)

  charValues = 'rgbcmywk'.';  %#'
  rgbValues = [eye(3); 1-eye(3); 1 1 1; 0 0 0];
  assert(~isempty(inColor),'convert_color:badInputSize',...
         'Input argument must not be empty.');

  if ischar(inColor)  %# Input is a character string

    [isColor,colorIndex] = ismember(inColor(:),charValues);
    assert(all(isColor),'convert_color:badInputContents',...
           'String input can only contain the characters ''rgbcmywk''.');
    outColor = rgbValues(colorIndex,:);

  elseif isnumeric(inColor) || islogical(inColor)  %# Input is a numeric or
                                                   %#   logical array
    assert(size(inColor,2) == 3,'convert_color:badInputSize',...
           'Numeric input must be an N-by-3 matrix');
    inColor = double(inColor);           %# Convert input to type double
    scaleIndex = max(inColor,[],2) > 1;  %# Find rows with values > 1
    inColor(scaleIndex,:) = inColor(scaleIndex,:)./255;  %# Scale by 255
    [isColor,colorIndex] = ismember(inColor,rgbValues,'rows');
    assert(all(isColor),'convert_color:badInputContents',...
           'RGB input must define one of the colors ''rgbcmywk''.');
    outColor = charValues(colorIndex(:));

  else  %# Input is an invalid type

    error('convert_color:badInputType',...
          'Input must be a character or numeric array.');

  end