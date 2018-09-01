function [XP,X,info,M] = mmatch_QP_PG(W,C,nFeature,k,varargin)

% INPUT -  W: m*m, pairwise matching result
%          C: 2*m, coordinate matrix
%          nFeature: number of feature candidates in each image
%          k: number of selected features
% OUTPUT - XP: m*m, joint matching result
%          X: m*k, mapping from image features to selected feature space
%          M: coordinates of selected features

tol_Y = 1e-4;
lambda = 1;% weight of geometric constraint
maxIter = 500;
maxIter_Y = 500;
verbose = true;

ivarargin = 1;
while ivarargin <= length(varargin)
    switch lower(varargin{ivarargin})
        case 'lambda'
            ivarargin = ivarargin+1;
            lambda = varargin{ivarargin};
        case 'maxiter'
            ivarargin = ivarargin+1;
            maxIter = varargin{ivarargin};
        case 'verbose'
            ivarargin = ivarargin+1;
            verbose = varargin{ivarargin};
        case 'rank'
            ivarargin = ivarargin+1;
            rank = varargin{ivarargin};
        otherwise
            fprintf('Unknown option ''%s'' is ignored!',varargin{ivarargin});
    end
    ivarargin = ivarargin+1;
end

m = size(W,1);
dimGroup = cumsum(nFeature);
dimGroup = [0; dimGroup(:)];
% randomly initialize Y
Y = rand(m,k);
t0 = tic;
% initialize Y by projected gradient descent
Y = initial_Y(Y,W,dimGroup);
% initialize X
U = zeros(m,k);
C_norm = C/((std(C(1,:))+std(C(2,:)))/2);
X = get_newX(U,Y,C_norm,0,k,rank,dimGroup,lambda);
%% update Y X Z
Rho = [1e0 1e1 1e2];
Iter = 0;
for i = 1:length(Rho)
    rho = Rho(i)
    for iter = 1:maxIter
        % update Y
        for iter_Y = 1:maxIter_Y
            tic;
            Y0 = Y;
            g = - W*Y + Y*(Y'*Y) + Rho(i)*(Y - X);
            st =  3*normest(Y'*Y,1e-2) + normest(W,1e-2) + Rho(i);
            Y = Y - g/st;
            for i2 = 1:length(dimGroup)-1
                ind = dimGroup(i2)+1:dimGroup(i2+1);
                Y(ind,:) = proj2dpam(Y(ind,:),1e-2);
            end
            t1 = toc;
            RelChg = norm(Y(:)-Y0(:))/sqrt(m);
            if verbose
                fprintf('Iter = %d, iter_Y = %d, Res = %d,t=(%.2f)\n',iter, iter_Y, RelChg, t1);
            end
            
            if  RelChg < tol_Y
                break
            end
        end
        % update X
        X0 = X;
        [X,M] = get_newX(Y,X,C_norm,Rho(i),k,rank,dimGroup,lambda);
        if sum(abs(X(:)-X0(:))) < 1
            break;
        end
    end
    Iter = Iter + iter;
end

%% overall match
info.time = toc(t0);
info.iter = Iter;
XP = X*X';
fprintf('Time = %fs #OverallIter = %d\n',info.time,info.iter);
end

function [X,M] = get_newX(Y,X,C,rho,K,rank,dimGroup,lambda)
%% geometric constraint
n = length(dimGroup) - 1;
for i = 1:n
    ind = dimGroup(i)+1:dimGroup(i+1);
    M(2*i-1:2*i,:) = C(:,ind)*X(ind,:);
end
%update Z
[U,S,V] = svd(M,'econ');
Z = U*S(:,1:rank)*V(:,1:rank)';

%update X
for i = 1:n
    ind = dimGroup(i)+1:dimGroup(i+1);
    Ci = C(:,ind);
    Zi = Z(2*i-1:2*i,:);
    Yi = Y(ind,:);
    % hungarian
    D = lambda*pdist2(Ci',Zi').^2;
    distMatrix =  D - 2*rho*Yi - min(min(D - 2*rho*Yi));
    assignment = assignmentoptimal(double(distMatrix));
    Xhi = zeros(length(ind),K);
    q = find(assignment >= 1);
    indices = sub2ind(size(Xhi), q, assignment(q));
    Xhi(indices) = 1;
    X(ind,:) = Xhi;
end
end

function Y = initial_Y(Y0, W, dimGroup)
fprintf('initialize Y...\n')
tol = 1e-4;
maxIter = 500;
lambda = [1e0 1e1 1e2 1e3 1e4];
Y = Y0;
m = size(W,1);
tic;
for i = 1:length(lambda)
    for iter = 1:maxIter
        Y0 = Y;
        % compute gradient and stepsize
        normr = 0;
        for j=1:length(dimGroup)-1
            ind = dimGroup(j)+1:dimGroup(j+1);
            Yi = Y(ind,:);
            YitYi = Yi'*Yi;
            Y_YYtd(ind,:) = Yi*YitYi;
            % stepsize of regularizer
            Hr = lambda(i)*(3*normest(Y(ind,:)'*Y(ind,:),1e-2) + 1);
            normr = max(Hr,normr);
        end
        
        gi = lambda(i)*(Y_YYtd - Y);
        gij = - W*Y + Y*(Y'*Y);
        g = gi + gij; %gradient
        st = 3*normest(Y'*Y) + normest(W,1e-2) + normr; %stepsize
        
        % update and project
        Y = Y - g/st;
        for k = 1:length(dimGroup)-1
            ind = dimGroup(k)+1:dimGroup(k+1);
            Y(ind,:) = proj2dpam(Y(ind,:),1e-2);
        end
        
        RelChg = norm(Y(:)-Y0(:))/sqrt(m);
        fprintf('lambda = %d, iter = %d, Res = %d\n',lambda(i), iter, RelChg);
        
        if  RelChg < tol
            break
        end
    end
end
t0 = toc;
fprintf('Y initialization finished, using %fs\n',t0);
end

