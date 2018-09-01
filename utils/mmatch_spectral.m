function [X,V,info] = mmatch_spectral(W,dimGroup,k)

k = min(k,size(W,1));
t0 = tic;
[V,~] = eigs(W,k,'la');
%  [V,~] = SVDS(W,'econ');
Y = rounding(V(:,1:k),dimGroup,0.5);
X = single(Y)*single(Y)'>0;
info.time = toc(t0);
                        
                        
                

    
    


