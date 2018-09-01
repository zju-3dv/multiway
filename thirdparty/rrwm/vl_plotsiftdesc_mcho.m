function [h h_b] =vl_plotsiftdesc_mcho(d,f,varargin)
% VL_PLOTSIFTDESCRIPTOR   Plot SIFT descriptor
%   VL_PLOTSIFTDESCRIPTOR(D) plots the SIFT descriptors D, stored as
%   columns of the matrix D. D has the same format used by VL_SIFT().
%
%   VL_PLOTSIFTDESCRIPTOR(D,F) plots the SIFT descriptors warped to
%   the SIFT frames F, specified as columns of the matrix F. F has the
%   same format used by VL_SIFT().
%
%   H=VL_PLOTSIFTDESCRIPTOR(...) returns the handle H to the line drawing
%   representing the descriptors.
%
%   REMARK. By default, the function assumes descriptors with 4x4
%   spatial bins and 8 orientation bins (Lowe's default.)
%
%   The function supports the following options
%
%   NumSpatialBins:: [4]
%     Number of spatial bins in each spatial direction.
%
%   NumOrientBins:: [8]
%     Number of orientation bis.
%
%   Magnif:: [3]
%     Magnification factor.
%
%   See also: VL_SIFT(), VL_PLOTFRAME(), VL_HELP().

% Copyright (C) 2007-12 Andrea Vedaldi and Brian Fulkerson.
% All rights reserved.
%
% This file is part of the VLFeat library and is made available under
% the terms of the BSD license (see the COPYING file).

% modified by Minsu Cho

magnif = 3.0 ;
NBO    = 8 ;
NBP    = 4 ;
maxv   = 0 ;
minv   = 0 ;
drawingtype = 'gradient';
bshowbox = true ; 
bContrastInsenstive = false;

if nargin > 1
  if ~ isnumeric(f)
    error('F must be a numeric type (use [] to leave it unspecified)') ;
  end
end

for k=1:2:length(varargin)
  opt=lower(varargin{k}) ;
  arg=varargin{k+1} ;
  switch opt
    case 'nbp'
      NBP = arg ;
    case 'nbo'
      NBO = arg ;
    case 'magnif'
      magnif = arg ;
    case 'maxv'
      maxv = arg ;
    case 'minv'
      minv = arg ;
    case 'drawingtype'
      drawingtype = arg;
    case 'showbox'
      bshowbox = arg ;  
    otherwise
      error(sprintf('Unknown option ''%s''', opt)) ;
  end
end

% --------------------------------------------------------------------
%                                                  Check the arguments
% --------------------------------------------------------------------

if size(d,1)*2 == NBP*NBP*NBO
    % contrast insenstive version
    NBO = NBO / 2;
    bContrastInsenstive = true;
end

if(size(d,1) ~= NBP*NBP*NBO)
  error('The number of rows of D does not match the geometry of the descriptor') ;
end

if nargin > 1
  if (~isempty(f) & size(f,1) ~= 6)
    error('F should be a 6xK matrix or the empty matrix');
  end

  if(~isempty(f) & size(f,2) ~= size(d,2))
    error('D and F have incompatible dimension') ;
  end
end

% Descriptors are often non-double numeric arrays
d = double(d) ;
K = size(d,2) ;

if nargin < 2 | isempty(f)
  f = repmat([0;0;1;0],1,K) ;
end

% --------------------------------------------------------------------
%                                                           Do the job
% --------------------------------------------------------------------

xg_all=[] ; yg_all=[] ;
xe_all=[] ; ye_all=[] ;
xb_all=[] ; yb_all=[] ;
h = [];   h_b = [];

for k=1:K
  %SBP = magnif * f(3,k) ;
  SBP = magnif; 
  %th=f(4,k) ;
  %c=cos(th) ;
  %s=sin(th) ;
  
  [xg, yg, xn, yn, xb, yb] = render_descr(d(:,k), NBP, NBO, maxv, minv, bContrastInsenstive) ;
  %xall = [xall SBP*(a*x-s*y)+f(1,k)] ;
  %yall = [yall SBP*(s*x+c*y)+f(2,k)] ;
  xg_all = [xg_all SBP*(f(3,k)*xg+f(5,k)*yg)+f(1,k)] ;
  yg_all = [yg_all SBP*(f(4,k)*xg+f(6,k)*yg)+f(2,k)] ;
  
  xe_all = [xe_all SBP*(f(3,k)*xn+f(5,k)*yn)+f(1,k)] ;
  ye_all = [ye_all SBP*(f(4,k)*xn+f(6,k)*yn)+f(2,k)] ;
  
  xb_all = [xb_all SBP*(f(3,k)*xb+f(5,k)*yb)+f(1,k)] ;
  yb_all = [yb_all SBP*(f(4,k)*xb+f(6,k)*yb)+f(2,k)] ;
end

if bshowbox
    h_b=line(xb_all,yb_all);
end
switch drawingtype
    case 'gradient'
        h=line(xg_all,yg_all);
    case 'edge'
        h=line(xe_all,ye_all);
    otherwise
end


% --------------------------------------------------------------------
function [xg,yg, xe, ye, xb,yb] = render_descr(d, BP, BO, maxv, minv, bContrastInsenstive)
% --------------------------------------------------------------------

[x,y] = meshgrid(-BP/2:BP/2,-BP/2:BP/2) ;

% Rescale d so that the biggest peak fits inside the bin diagram
if maxv
    %d = 0.5 * (d - minv) / (maxv - minv +eps);
    d = 0.5 * d / (maxv + eps) ;
else
    %d = 0.5 * (d - min(d(:))) / (max(d(:)) - min(d(:)) +eps ) ;
    d = 0.5 * d / (max(abs(d(:))) +eps) ;
end

% supress the negative values
d( d < 0) = 0;

% We have BP*BP bins to plot. Here are the centers:
xc = x(1:end-1,1:end-1) + 0.5 ;
yc = y(1:end-1,1:end-1) + 0.5 ;

% We scramble the the centers to have the in row major order
% (descriptor convention).
xc = xc' ;
yc = yc' ;

% Each spatial bin contains a star with BO tips
xc = repmat(xc(:)',BO,1) ;
yc = repmat(yc(:)',BO,1) ;

% Do the stars
if ~bContrastInsenstive
    th=linspace(0,2*pi,BO+1) ;
else
    th=linspace(0,pi,BO+1) ;
end
th=th(1:end-1) ;
xd = repmat(cos(th), 1, BP*BP) ;
yd = repmat(sin(th), 1, BP*BP) ;
xd = xd .* d(:)' ;
yd = yd .* d(:)' ;

% Re-arrange in sequential order the lines to draw
nans = NaN * ones(1,BP^2*BO) ;
x1 = xc(:)' ;
y1 = yc(:)' ;
x2 = xc(:)' + xd ;
y2 = yc(:)' + yd ;
xstars = [x1;x2;nans] ;
ystars = [y1;y2;nans] ;

% Do the edge features normal to the stars 
xde = repmat(cos(th+pi/2), 1, BP*BP) ;
yde = repmat(sin(th+pi/2), 1, BP*BP) ;
xde = xde .* d(:)' ;
yde = yde .* d(:)' ;

% Re-arrange in sequential order the lines to draw
x1e = xc(:)' - xde;
y1e = yc(:)' - yde;
x2e = xc(:)' + xde ;
y2e = yc(:)' + yde ;
xstarse = [x1e;x2e;nans] ;
ystarse = [y1e;y2e;nans] ;

% Horizontal lines of the grid
%nans = NaN * ones(1,BP+1);
%xh = [x(:,1)' ; x(:,end)' ; nans] ;
%yh = [y(:,1)' ; y(:,end)' ; nans] ;
xh = [x([1 end],1)' ; x([1 end],end)' ; NaN NaN] ;
yh = [y([1 end],1)' ; y([1 end],end)' ; NaN NaN] ;

% Verical lines of the grid
%xv = [x(1,:) ; x(end,:) ; nans] ;
%yv = [y(1,:) ; y(end,:) ; nans] ;
xv = [x(1,[1 end]); x(end,[1 end]); NaN NaN] ;
yv = [y(1,[1 end]); y(end,[1 end]); NaN NaN] ;

xg=[xstars(:)'] ;
yg=[ystars(:)'] ;

xe=[xstarse(:)'] ;
ye=[ystarse(:)'] ;

xb=[xh(:)' xv(:)'] ;
yb=[yh(:)' yv(:)'] ;
