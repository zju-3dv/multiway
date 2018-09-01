function drawCircle(x, y, r, S, varargin)

linewidth = 1;
for k=1:2:length(varargin)
  opt=lower(varargin{k}) ;
  arg=varargin{k+1} ;
  switch opt
    case 'linewidth'
      linewidth = arg ;
    otherwise
      error(sprintf('Unknown option ''%s''', opt)) ;
  end
end

theta = 0 : (2 * pi / 100) : (2 * pi);
pline_x = r * cos(theta) + x;
pline_y = r * sin(theta) + y;

plot(pline_x, pline_y, S, 'linewidth', linewidth);