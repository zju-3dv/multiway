% Least Squares Affine Transformation
% ELEC 301 Group Project
% 11/29/2009
% Jeffrey Bridge, Robert Brockman II, Stamatios Mastrogiannis
%
% Calculate the least squares affine transformation for two corresponding
% sets of pixel locations.
% px inputs are of the form:
%[ x_1 y_1
%  x_2 y_2
%  :   :
%  x_N y_N ]
%
% [x'] = [a, b] * [x] + [e]
% [y']   [c, d]   [y]   [f]

function Aff = l2aff(pxold, pxnew)
    b = reshape(pxnew.', [], 1);
    A = makenice(pxold);
    x = pinv(A) * b; % Was psinv, our version of computing the pseudoinv
    Aff = [x(1), x(2), x(5); ...
          x(3), x(4), x(6)];
return

function A = makenice(pxold)
    [r, c] = size(pxold);
    A = zeros(2*r, 6);
    for k=1:r
        x = pxold(k,1);
        y = pxold(k,2);
        %correspond to a, b, c, d, e, f
        A(2*k-1, :) = [x, y, 0, 0, 1, 0];
        A(2*k  , :) = [0, 0, x, y, 0, 1];
    end
return