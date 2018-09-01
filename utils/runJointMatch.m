function [jMatch,jmInfo,timeInfo] = runJointMatch(pMatch,C,varargin)

% INPUT -  pMatch:  pairwise matching result
%          C:       2*m coordinate matrix,
%                   where m is the total number of keypoints in the image collection
% OUTPUT - jMatch:  joint matching result

% parameter settings
Size = []; % number of selected features
pSelect = 1; % ratio of selection
wSparsity = 0; % weight of sparsity
method = 'pg'; % joint matching method
rank = 3; % rank for geometric constraint
lambda = 1; % the weight of geometric constraint

ivargin = 1;
while ivargin <= length(varargin)
    switch lower(varargin{ivargin})
        case 'method'
            ivargin = ivargin + 1;
            method = varargin{ivargin};
        case 'univsize'
            ivargin = ivargin + 1;
            Size = varargin{ivargin};
        case 'pselect'
            ivargin = ivargin + 1;
            pSelect = varargin{ivargin};
        case 'wsparsity'
            ivargin = ivargin + 1;
            wSparsity = varargin{ivargin};
        case 'lambda'
            ivargin = ivargin + 1;
            lambda = varargin{ivargin};
        case 'rank'
            ivargin = ivargin + 1;
            rank = varargin{ivargin};
        otherwise
            fprintf('Unknown option ''%s'' is ignored !!!\n',varargin{ivargin});
    end
    ivargin = ivargin + 1;
end

nFeature = zeros(size(pMatch,1),1);
filename = cell(size(pMatch,1),1);
nMatches = 0;
for i = 1:size(pMatch,1)
    for j = i+1:size(pMatch,2)
        if ~isempty(pMatch(i,j).X)
            nMatches = nMatches + sum(pMatch(i,j).X);
        end
        if ~isempty(pMatch(i,j).nFeature)
            nFeature(i) = pMatch(i,j).nFeature(1);
            nFeature(j) = pMatch(i,j).nFeature(2);
        end
        if ~isempty(pMatch(i,j).filename)
            filename(i) = pMatch(i,j).filename(1);
            filename(j) = pMatch(i,j).filename(2);
        end
    end
end

cumIndex = cumsum([0;nFeature]);
ind1 = zeros(1,nMatches);
ind2 = zeros(1,nMatches);
flag = zeros(1,nMatches);
score = zeros(1,nMatches);
z = 0;
for i = 1:size(pMatch,1)
    for j = i+1:size(pMatch,2)
        try
            matchList = double(pMatch(i,j).matchInfo.match);
            n = length(pMatch(i,j).X);
            ind1(z+1:z+n) = cumIndex(i) + matchList(1,:);
            ind2(z+1:z+n) = cumIndex(j) + matchList(2,:);
            flag(z+1:z+n) = pMatch(i,j).X;
            score(z+1:z+n) = mat2gray(pMatch(i,j).Xraw);
            z = z + n;
        catch
        end
    end
end

% original scores
M = sparse(ind1,ind2,score,cumIndex(end),cumIndex(end));
M = M + M';
% binary scores
Mbin = sparse(ind1,ind2,flag,cumIndex(end),cumIndex(end));
Mbin = Mbin + Mbin';
vM = Mbin;

Size = min(Size,min(nFeature));
Z = [];
fprintf('>> Run joint match, problem size = (%d,%d) ...\n',size(vM,1),Size);
switch lower(method)
    case 'spectral'
        % Spectral method
        [M_out,eigV,timeInfo] = mmatch_spectral(vM,nFeature,Size);
    case 'matchlift'
        % MatchLift
        [M_out,eigV,timeInfo] = mmatch_MatchLift(vM,nFeature,Size);
    case 'als'
        % MatchALS
        [M_out,eigV,timeInfo] = mmatch_CVX_ALS(vM,nFeature,'univsize',Size,'pSelect',pSelect,...
            'tol',5e-4,'beta',wSparsity);
    case 'pg'
        % the proposed method
        [M_out,eigV,timeInfo,Z] = mmatch_QP_PG(vM,C,nFeature,Size,'lambda',lambda,'rank',rank);
end

csum = cumsum([0;nFeature]);
for i = 1:size(pMatch,1)
    for j = i+1:size(pMatch,2)
        [ind1,ind2] = find(M_out(csum(i)+1:csum(i+1),csum(j)+1:csum(j+1)));
        if ~isempty(ind1)
            % remove conflict
            Xraw = vM(csum(i)+1:csum(i+1),csum(j)+1:csum(j+1));
            Xraw = Xraw(sub2ind(size(Xraw),ind1,ind2));
            X = greedyMatch([ind1';ind2'],Xraw);
            % store results
            jMatch(i,j).matchInfo.match = [ind1';ind2'];
            jMatch(i,j).X = X;
            jMatch(i,j).Xraw = Xraw;
            jMatch(i,j).filename = [filename(i),filename(j)];
            jMatch(i,j).nFeature = [nFeature(i),nFeature(j)];
        end
    end
end

jmInfo.eigV = eigV;
jmInfo.nFeature = csum;
jmInfo.filename = filename;
jmInfo.time = timeInfo.time;
jmInfo.Z = Z;