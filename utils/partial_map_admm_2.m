% ADMM Algorithm
% Para example:
% Para.nIterations = 1000;
% Para.mu_init = 0.1;
% Para.rho = 1.008;

function [X,info] = partial_map_admm_2(OPT)

Para.nIterations = 1000; 
Para.mu_init = 0.1; 
Para.rho = 1.008;

if isfield(OPT,'maxIter')
    Para.nIterations = OPT.maxIter;
end

dim = size(OPT.W, 1);

C = ones(dim, dim) - 2*OPT.W;

% Equality constraints
A = OPT.A;
b = OPT.b;
nr = length(b);
%A = A(2:nr, :);
%b = b(2:nr);
%A = A(1:(2*dim-1), :);
%b = b(1:(2*dim-1), :);


% Allocate variables

X = zeros(dim, dim);
S = X;
Z = X;
y1 = zeros(length(b), 1);

mu = Para.mu_init;

t0 = tic;
for iteration = 1:Para.nIterations
    % Optimizing y1
    TP = -C + S + mu*X + Z;
    rhs = A*reshape(TP, [dim*dim,1]) - mu*b;
    y1 = (A*A')\rhs;
    
    
    % Optimizing Z
    TP = reshape(A'*y1, [dim,dim]);
    TP = TP + TP' -diag(diag(TP));
    TP = C + TP - mu*X;
    Z = max(TP - S, 0);
    
    % Optimizing S
    TP = TP - Z;
    TP = (TP + TP')/2;
    try
        [U,V] = eig(TP);
    catch 
        [U,V] = laneig(TP);
    end
    ids = find(diag(V) > 0);
    S = U(:, ids)*V(ids, ids)*U(:, ids)';
    
    % Updating X
    X = (S-TP)/mu;
    X = (X+X')/2;
    
    mu = mu*Para.rho;
    if mod(iteration, 50) == 0
        res1 = norm(A*reshape(X, [dim*dim,1]) - b);
        e = norm(res1);
%         fprintf('   res1 = %f, mu = %f, e= %f\n', res1, mu, e);
        if e < 1e-3
            break;
        end
    end
end

info.time = toc(t0);
info.iter = iteration;