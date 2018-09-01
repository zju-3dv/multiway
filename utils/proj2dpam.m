function X = proj2dpam(Y,tol)

if nargin < 2
    tol = 1e-4;
end

X0 = Y;
X = Y;
I2 = 0;

for iter = 1:10

    % row projection
    X1 = projR(X0+I2); I1 = X1-(X0+I2);
    X2 = projC(X0+I1); I2 = X2-(X0+I1);
    
    chg = sum(abs(X2(:)-X(:)))/numel(X);
%     fprintf('%f\n',chg);
    X = X2;
    if chg < tol
        return
    end

end

end

function X = projR(X)
    for i = 1:size(X,1)
        X(i,:) = proj2pav(X(i,:)')';
    end
end

function X = projC(X)
    for j = 1:size(X,2)
        X(:,j) = proj2pavC(X(:,j));
    end
end

function x = proj2pav(y)

y(y<0) = 0;
if sum(y) < 1
    x = y;
else
    u = sort(y,'descend');
    sv = cumsum(u);
    rho = find(u >(sv-1)./(1:length(u))',1,'last');
    theta = max(0,(sv(rho)-1)/rho);
    x = max(y-theta,0);
end
end
function x = proj2pavC(y)
% project an n-dim vector y to the simplex Dn
% Dn = { x : x n-dim, 1 >= x >= 0, sum(x) = 1}

% (c) Xiaojing Ye
% xyex19@gmail.com
%
% Algorithm is explained as in the linked document
% http://arxiv.org/abs/1101.6081
% or
% http://ufdc.ufl.edu/IR00000353/
%
% Jan. 14, 2011.

m = length(y); bget = false;

s = sort(y,'descend'); tmpsum = 0;

for ii = 1:m-1
    tmpsum = tmpsum + s(ii);
    tmax = (tmpsum - 1)/ii;
    if tmax >= s(ii+1)
        bget = true;
        break;
    end
end
    
if ~bget, tmax = (tmpsum + s(m) -1)/m; end;

x = max(y-tmax,0);

end