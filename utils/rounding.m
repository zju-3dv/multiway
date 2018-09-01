function Y = rounding(A,dimGroup,th)

% Input:
% Y - eigen map of points
% dimGroup: number of points in each group
% Output: rounded assignment matrix

if nargin < 3
    th = 0.5;
end

% normalize in order to calculate correlation coefficient of points
A = normr(A);

N = cumsum(dimGroup);
flag = false(size(A,1),1); % indicator if a point has been assigned
Y = false(0,0);

for i = 1:size(A,1)
    if flag(i) == false
        p = size(Y,2)+1;
        Y(i,p) = 1;
        flag(i) = true;
        ob = find(N>=i,1,'first');
        if ob < length(N)
            for ob1 = ob+1:length(N)
                z = [];
                cc = A(i,:)*A(N(ob1-1)+1:N(ob1),:)';
                mx = max(cc(flag(N(ob1-1)+1:N(ob1))==false));
                if mx > th
                    z = find(cc==mx) + N(ob1-1);
                end
                if ~isempty(z)
                    Y(z,p) = 1;
                    flag(z) = true;
                end
            end
        end
    end
end