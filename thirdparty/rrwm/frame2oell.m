% --------------------------------------------------------------------
function eframes = frame2oell(frames)
% FRAMES2OELL  Convert generic frame to oriented ellipse
%   EFRAMES = FRAME2OELL(FRAMES) converts the frames FRAMES to
%   oriented ellipses EFRAMES. This is useful because many tasks are
%   almost equivalent for all kind of regions and are immediately
%   reduced to the most general case.

% Determine the kind of frames
[D,K] = size(frames) ;
switch D
  case 2, kind = 'point' ;
  case 3, kind = 'disk' ;
  case 4, kind = 'odisk' ;
  case 5, kind = 'ellipse' ;
  case 6, kind = 'oellipse' ;
  otherwise
    error(['FRAMES format is unknown']) ;
end
eframes = zeros(6,K) ;

% Convert frames to oriented ellipses
switch kind
  case 'point'
    eframes(1:2,:) = frames(1:2,:) ;

  case 'disk'
    eframes(1:2,:) = frames(1:2,:) ;
    eframes(3,:)   = frames(3,:) ;
    eframes(6,:)   = frames(3,:) ;

  case 'odisk'
    r = frames(3,:) ;
    c = r.*cos(frames(4,:)) ;
    s = r.*sin(frames(4,:)) ;

    eframes(1:2,:) = frames(1:2,:) ;
    eframes(3:6,:) = [c ; s ; -s ; c] ;

  case 'ellipse'
    eframes(1:2,:) = frames(1:2,:) ;
    eframes(3:6,:) = mapFromS(frames(3:5,:)) ;

  case 'oellipse'
    eframes = frames ;
end

% --------------------------------------------------------------------
function A = mapFromS(S)
% --------------------------------------------------------------------
% Returns the (stacking of the) 2x2 matrix A that maps the unit circle
% into the ellipses satisfying the equation x' inv(S) x = 1. Here S
% is a stacked covariance matrix, with elements S11, S12 and S22.

tmp = sqrt(S(3,:)) + eps ;
A(1,:) = sqrt(S(1,:).*S(3,:) - S(2,:).^2) ./ tmp ;
A(2,:) = zeros(1,length(tmp));
A(3,:) = S(2,:) ./ tmp ;
A(4,:) = tmp ;