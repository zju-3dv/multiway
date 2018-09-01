function [Rsh,Rp] = initRotation(W,K)

wm = mean(W,2);
W = W - wm*ones(1,size(W,2));
F = size(W,1)/2;
P = size(W,2);

%% SVD
if (2*F)>P
    [U,D,V] = svd(W,0);
else
    [V,D,U] = svd(W',0);
end;

Mhat = U(:,1:3*K)*sqrt(D(1:3*K,1:3*K));
K3 = 3*K;%P; %% CHECK THISS
MhatP = U(:,1:K3)*sqrt(D(1:K3,1:K3));
%%  Compute rotations

PI_Hat = Mhat;
nPose = F;
A = zeros(2*nPose,(3*K)*(3*K));
A_Hat = zeros(2*nPose,(3*K)*(3*K+1)/2);

for f=1:nPose
    PI_Hat_f = PI_Hat(2*f-1:2*f,:);
    AA = kron(PI_Hat_f,PI_Hat_f);
    
    A(2*f-1,:) = AA(1,:) - AA(4,:);
    A(2*f-0,:) = AA(2,:);
    
    count = 0;
    for i=1:3*K
        for j=i:3*K
            count = count+1;
            if(i==j)
                A_Hat(2*f-1,count)=A(2*f-1,(3*K)*(i-1)+j);
                A_Hat(2*f-0,count)=A(2*f-0,(3*K)*(i-1)+j);
            else
                A_Hat(2*f-1,count)=A(2*f-1,(3*K)*(i-1)+j)+A(2*f-1,(3*K)*(j-1)+i);
                A_Hat(2*f-0,count)=A(2*f-0,(3*K)*(i-1)+j)+A(2*f-0,(3*K)*(j-1)+i);
            end
        end
    end
end

[U_A D_A V_A] = svd(A_Hat);

for i = 1:size(V_A,2)
    if(V_A(1,i) < 0)
        V_A(:,i) = -V_A(:,i);
    end
end

X_Line = V_A(:,count-(2*K^2-K)+1:end); %This is the null space

b = ones(2*K^2-K,1);
cvx_begin quiet
variable t(2*K^2-K,1)
variable Q_Hat(3*K,3*K) symmetric
variable XX(3*K,3*K) symmetric
variable YY(3*K,3*K) symmetric
Q_Line = X_Line*t;
Z = [XX Q_Hat;
    Q_Hat' YY];      % the constructed matrix
minimize ( 1/2*(trace(XX) + trace(YY)))
subject to
Z==semidefinite(6*K,6*K);
b'*t==1;
Q_Line = X_Line*t;
count = 0;
for i = 1:3*K
    for j = i:3*K
        count = count + 1;
        Q_Hat(i,j) == Q_Line(count);
        Q_Hat(j,i) == Q_Line(count);
    end
end
cvx_end

scale = b'*t;
Q_Hat = double(Q_Hat/scale);
[U_Hat D_Hat V_Hat] = svd(Q_Hat);
G_Hat = U_Hat(:,1:3)*sqrt(D_Hat(1:3,1:3));
%Nonlinear optimization
options = optimset('Display','off','MaxFunEval',200000,'MaxIter',1000,'TolFun',1e-4,'TolX',1e-4);
warning off;
[G, fval] = fminunc(@evalQ_regularization,G_Hat,options,PI_Hat);  %fminunc  fminsearch
[Rsh Rp] = Recover_Rotation(PI_Hat,G);

end


function [R_Recover Rsh] = Recover_Rotation(PI_Hat,G)

nPose = size(PI_Hat,1)/2;  %Number of frames

R_Recover = zeros(2*nPose,3); %Recovered matrix 2F*3

c_fk1 = zeros(nPose,1);
c_fk2 = zeros(nPose,1);

for f=1:nPose
    Eq = PI_Hat(2*f-1:2*f,:)*G*G'*PI_Hat(2*f-1:2*f,:)';

    c_fk1(f) = sqrt(abs(Eq(1,1))); %be determined by the rotation
    c_fk2(f) = sqrt(abs(Eq(2,2))); %be determined by the rotation
    
    R_Recover(2*f-1,:) = PI_Hat(2*f-1,:)*G/(c_fk1(f));
    R_Recover(2*f-0,:) = PI_Hat(2*f-0,:)*G/(c_fk2(f));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Here we have obtained the rotation for each frame and the coefficient,
%however there is a sign ambiguity, then the problem is how to fix it, the
%constraint between two frames is that the angle is less than 90 degrees.
r1 = R_Recover(1,:);
r2 = R_Recover(2,:);
r3 = cross(r1,r2);
R_f_1 = [r1;r2;r3];

if(det(R_f_1) < 0)
    R_f_1 = - R_f_1;
end

R_Recover(1:2,:) = R_f_1(1:2,:);

theta1 = zeros(nPose,1);
theta2 = zeros(nPose,1);

for f = 2:nPose
    r1 = R_Recover(2*f-1,:);
    r2 = R_Recover(2*f-0,:);
    r3 = cross(r1,r2);
    R_f = [r1;r2;r3];
    
    if(det(R_f) < 0)
        R_f = - R_f;
    end
        
    theta1(f) = acos((trace(R_f'*R_f_1)-1)/2)*180/pi;
    
    if(theta1(f) > 90)
        R_f(1:2,:) = - R_f(1:2,:);
    end
    
    theta2(f) = acos((trace(R_f'*R_f_1)-1)/2)*180/pi;
        
    R_Recover(2*f-1:2*f,:) = R_f(1:2,:);
    
    R_f_1 = R_f;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rsh is the big R 2F*3F diagonal
for i = 1:nPose
    Rsh(2*(i-1)+1:2*(i-1)+2,3*(i-1)+1:3*(i-1)+3) = R_Recover(2*i-1:2*i,:);
end
end
