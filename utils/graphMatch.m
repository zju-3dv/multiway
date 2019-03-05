function mdata = graphMatch(viewPair,varargin)

methodGM = 'rw'; % method for solving graph matching
methodDisc = 'greedy'; % method for discretization 
wEdge = 1; % weigth of rigidicity 
wAngle = 0.3;
kNNInit = 50; % number of candidate correpondences for each feature
nMaxInit = 5000; % total number of candidate correpondences after optimization
thScore = 0; % threshold of matching scores below which the match is ignored
thRatio = 1; 
thDist = 0;

ivargin = 1;
while ivargin <= length(varargin)
switch(lower(varargin{ivargin}))
    case 'methodgm'
        ivargin = ivargin+1;
        methodGM = varargin{ivargin};
    case 'methoddisc'
        ivargin = ivargin+1;
        methodDisc = varargin{ivargin};
    case 'wedge'
        ivargin = ivargin+1;
        wEdge = varargin{ivargin};
    case 'knninit'
        ivargin = ivargin+1;
        kNNInit = varargin{ivargin};
    case 'nmaxinit'
        ivargin = ivargin+1;
        nMaxInit = varargin{ivargin};
    case 'thscore'
        ivargin = ivargin+1;
        thScore = varargin{ivargin};
    case 'thratio'
        ivargin = ivargin+1;
        thRatio = varargin{ivargin};
    case 'thdist'
        ivargin = ivargin+1;
        thDist = varargin{ivargin};
    case 'paramscript'
        ivargin = ivargin + 1;
        paramScript = varargin{ivargin};
        eval(paramScript);
    otherwise
        fprintf('Unknown option ''%s'' is ignored!\n',varargin{ivargin});
end
ivargin = ivargin+1;
end

% obtain initial matches
fprintf('- looking for initial matches ...\n');
matchInfo = initMatch(viewPair,kNNInit,nMaxInit,thScore,thRatio,thDist);
match = matchInfo.match;
mdata.matchInfo = matchInfo;
mdata.nFeature = [viewPair.nfeature];

if isempty(match)
    mdata.Xraw = [];
    mdata.X = [];
    return
end

% graph matching to introduce geometric constraint
if wEdge <= 0
    Xraw = double(matchInfo.sim); 
else
    % construct graphs
    fprintf('- constructing graphs ...\n');
    [uniq_feat1,~,new_feat1] = unique(match(1,:),'sorted');
    [uniq_feat2,~,new_feat2] = unique(match(2,:),'sorted');
    cand_matchlist_uniq = [new_feat1';new_feat2'];
    tmp = viewPair(1).frame(1:2,uniq_feat1);
    tmp = tmp / mean(std(tmp,1,2));
    tmp = [tmp;[1,0,0,1]'*ones(1,size(tmp,2))];
    edgeAttr1 = computeEdgeAttr(double(tmp));
    edgeAttr1(2,:) = edgeAttr1(2,:)*wAngle;
    tmp = viewPair(2).frame(1:2,uniq_feat2);
    tmp = tmp / mean(std(tmp,1,2));
    tmp = [tmp;[1,0,0,1]'*ones(1,size(tmp,2))];
    edgeAttr2 = computeEdgeAttr(double(tmp));
    edgeAttr2(2,:) = edgeAttr2(2,:)*wAngle;
    fE1 = makeFeatBin(edgeAttr1);
    fE2 = makeFeatBin(edgeAttr2);
    eSimVal = computeDotProdSimilarity_sym( int32(cand_matchlist_uniq), fE1, fE2); % symmetric affinities
    vSimVal = matchInfo.sim;
    mSimVal = repmat(vSimVal,[numel(vSimVal),1]); % distribute the v score
    affinityMatrix = wEdge*eSimVal + mSimVal + mSimVal';
    affinityMatrix(1:size(affinityMatrix,1)+1:end) = 0;
%     nlink = min(length(uniq_feat1),length(uniq_feat2));
%     affinityMatrix = wEdge*eSimVal/nlink + diag(vSimVal);
    affinityMatrix(affinityMatrix<1e-3) = 0;
    affinityMatrix = sparse(double(affinityMatrix));

    % choose graph matching method
    switch lower(methodGM)
        case 'sm'
            % spectral matching
            fprintf('- running spectral matching ...\n');
            [Xraw,~] = eigs(affinityMatrix,1,'LM'); 
            Xraw = abs(Xraw);
        case 'rw'
            % random walk
            fprintf('- running random walk matching ...\n');
            [group1,group2] = make_group12(match(1:2,:));
            Xraw = RRWM(affinityMatrix,group1,group2);
        otherwise
            fprintf('Unknown methodGM for graph matching!\n');
            Xraw = double(matchInfo.sim);     
    end
end

% discretization
switch lower(methodDisc)
    case 'greedy'
        X = greedyMatch(match,Xraw);
    case 'hungary'
        X = optimalMatch(match,Xraw);
    otherwise
        sprintf('!!! Unknown methodDisc for discretization.\n');
        X = greedyMatch(match,Xraw);
end
try
    mdata.Xraw = single(Xraw);
    mdata.X = X;
catch
    pause
end
end
