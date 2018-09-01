function X = spectralMaching( A )

[X,~] = eigs(A,1,'la'); 
X = abs(X);

end

