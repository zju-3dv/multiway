function y = proj2kav(x,K)

opts.Display = 'off';
opts.TolX = 1e-6;
x = double(x);
y = quadprog(eye(length(x)),-x,[],[],ones(1,length(x)),K,...
        zeros(size(x)),ones(size(x)),[],opts);