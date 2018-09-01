function [X_bin,A,info] = mmatch_CVX_ALS(W,dimGroup,varargin)
% This function is to solve
% min - <W,X> + alpha||X||_* + beta||X||_1, st. X \in C

% The problem is rewritten as
% <beta-W,AB^T> + alpha/2||A||^2 + alpha/2||B||^2
% st AB^T=Z, Z\in\Omega

% ---- Output:
% X: a sparse binary matrix indicating correspondences
% A: AA^T = X;
% info: other info.

% ---- Required input:
% W: sparse input matrix storing scores of pairwise matches
% dimGroup: a vector storing the number of points on each objects

% ---- Other options:
% maxRank: the restricted rank of X* (select it as large as possible)
% alpha: the weight of nuclear norm
% beta: the weight of l1 norm
% pSelect: propotion of selected points, i.e., m'/m in section 5.4 in the paper
% tol: tolerance of convergence
% maxIter: maximal iteration
% verbose: display info or not
% eigenvalues: output eigenvalues or not

% optional paramters
alpha = 50;
beta = 0.1;
maxRank = max(dimGroup)*4;
pSelect = 1;
tol = 5e-4;
maxIter = 1000;
verbose = true;
eigenvalues = false;

ivarargin = 1;
while ivarargin <= length(varargin)
    switch lower(varargin{ivarargin})
        case 'dimgroup'
            ivarargin = ivarargin+1;
            dimGroup = varargin{ivarargin};
        case 'univsize'
            ivarargin = ivarargin+1;
            maxRank = varargin{ivarargin};
        case 'alpha'
            ivarargin = ivarargin+1;
            alpha = varargin{ivarargin};
        case 'beta'
            ivarargin = ivarargin+1;
            beta = varargin{ivarargin};
        case 'pselect'
            ivarargin = ivarargin+1;
            pSelect = lower(varargin{ivarargin});
        case 'tol'
            ivarargin = ivarargin+1;
            tol = varargin{ivarargin};
        case 'maxiter'
            ivarargin = ivarargin+1;
            maxIter = varargin{ivarargin};
        case 'verbose'
            ivarargin = ivarargin+1;
            verbose = varargin{ivarargin};
        case 'eigenvalues'
            ivarargin = ivarargin+1;
            eigenvalues = varargin{ivarargin};
        otherwise
            fprintf('Unknown option ''%s'' is ignored!',varargin{ivargin});
    end
    ivarargin = ivarargin+1;
end

fprintf('Running MatchALS: alpha = %.2f, beta = %.2f, maxRank = %d, pSelect = %.2f \n',...
    alpha,beta,maxRank,pSelect);

W(1:size(W,1)+1:end) = 0;
W = (W+W')/2;
X = single(full(W));
Z = single(full(W));
Y = single(zeros(size(X)));
mu = 64;

n = size(X,1);
maxRank = min(n,ceil(maxRank));
A = rand(n,maxRank);

dimGroup = cumsum(dimGroup);
dimGroup = [0;dimGroup(:)];

t0 = tic;
for iter = 1:maxIter
    
    X0 = X;
    X = Z - (double(Y)-W+beta)/mu;
    B = ((A'*A+alpha/mu*eye(maxRank))\(A'*X))';
    A = ((B'*B+alpha/mu*eye(maxRank))\(B'*X'))';
    X = A*B';
    
    Z = X + Y/mu;
    diagZ = diag(Z);
    % enforce the self-matching to be null
    for i = 1:length(dimGroup)-1
        ind1 = dimGroup(i)+1:dimGroup(i+1);
        Z(ind1,ind1) = zeros(length(ind1),length(ind1)); 
    end
    % optimize for diagnal elements
    if pSelect == 1
        Z(1:size(Z,1)+1:end) = 1;
    else
        diagZ = proj2kav(diagZ,pSelect*length(diagZ));
        Z(1:size(Z,1)+1:end) = diagZ;
    end
    % rounding all elements to [0,1]
    Z(Z<0) = 0;
    Z(Z>1) = 1;
    
    Y = Y + mu*(X-Z);
    
    pRes = norm(X(:)-Z(:))/n;
    dRes = mu*norm(X(:)-X0(:))/n;
    if verbose 
        fprintf('Iter = %d, Res = (%d,%d), mu = %d \n',iter,pRes,dRes,mu);
    end
    
    if  pRes < tol && dRes < tol
        break
    end
    
    if pRes>10*dRes
        mu = 2*mu;
    elseif dRes>10*pRes
        mu = mu/2;
    else
    end
    
end

X = (X+X')/2;

info.time = toc(t0);
info.iter = iter;
if eigenvalues
    info.eigenvalues = eig(X);
end

X_bin = sparse(X>0.5);

fprintf('Alg terminated. Time = %d, #Iter = %d, Res = (%d,%d), mu = %d \n',...
    info.time,info.iter,pRes,dRes,mu);
