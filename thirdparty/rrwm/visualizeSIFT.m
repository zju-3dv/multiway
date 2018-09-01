function visualizeSIFT(d, f, varargin)
% visualize SIFT feature like HOG ( negative values allowed )
% by Minsu Cho

descScale = 2.0;
NBO = 8;

NBP = sqrt(size(d,1)/NBO);
if NBP*NBP*NBO ~= size(d,1) % contrast insensitive
    NBP = sqrt(2*size(d,1)/NBO);
end

if nargin > 1
  if ~ isnumeric(f)
    error('F must be a numeric type (use [] to leave it unspecified)') ;
  end
end

for k=1:2:length(varargin)
  opt=lower(varargin{k}) ;
  arg=varargin{k+1} ;
  switch opt
    case 'descscale'
      descScale = arg;  
      %magnif = arg*0.5;
    case 'offset'
      f(1,:) = f(1,:) + arg(1);
      f(2,:) = f(2,:) + arg(2);
    %case 'magnif'
    %  magnif = arg ;
    otherwise
      error(sprintf('Unknown option ''%s''', opt)) ;
  end
end

magnif = 0.5*descScale*4.0/NBP; % 0.5 * descScale *4.0/NBP (descScale)
%maxv = double(max(max(d))) ;
%minv = double(min(min(d))) ;

f = frame2oell(f);
%vl_plotframe([f(1:2,:); f(3:6,:)], 'color', 'y');
[h_grad h_box] = vl_plotsiftdesc_mcho(d, f, 'drawingtype', 'gradient', 'magnif', magnif, 'NBO', NBO, 'NBP', NBP) ;
[h_edge] = vl_plotsiftdesc_mcho(d, f, 'drawingtype', 'edge', 'magnif', magnif, 'NBO', NBO, 'NBP', NBP) ;
%h3 = vl_plotsiftdesc_grad(d, f, 'maxv', maxv, 'minv', minv) ;
set(h_grad,'color','g', 'linewidth',1) ;
set(h_box,'color','b', 'linewidth',1) ;
set(h_edge,'color','r', 'linewidth',2) ;