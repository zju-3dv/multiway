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

% if orientation info exists
if size(S,1) == 4
    
end