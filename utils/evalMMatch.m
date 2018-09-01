function [overlap,precision,recall] = evalMMatch(X,Xgt)

X = X > 0;
Xgt = Xgt > 0;

s1 = triu(X,1);
s2 = triu(Xgt,1);

overlap = nnz(s1&s2)/(nnz(s1|s2)+eps);
precision = nnz(s1&s2)/(nnz(s1)+eps);
recall = nnz(s1&s2)/(nnz(s2)+eps);