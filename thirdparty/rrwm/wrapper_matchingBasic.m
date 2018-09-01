function [ res ] = wrapper_matchingBasic( method, cdata )

aparam.bOneToOne = 1;
bShow = 0; 

if 0
    if bShow, hFig1 = figure(481);    clf; end
    set(0,'CurrentFigure',hFig1);   imshow(appendimages(cdata.view(1).img,cdata.view(2).img));
    t = clock; time_tag = sprintf('%04d_%02d_%02d_%02d_%02d_%02d',t(1), t(2), t(3), t(4), t(5), round(t(6)));
    saveas(hFig1, sprintf('./save/%s__input.png',time_tag));
    %showInitialMatches2; %saveas(hFig1, sprintf('./save/%s__initial.png',time_tag));
end

% command to call a matching function
strFunc = ['feval(@' func2str(method.fhandle)];
for j = 1:length(method.variable)
    strFunc = [ strFunc ', cdata.' method.variable{j} ];
end
strFunc = [strFunc ')'];

res.score_GM = -9999;   res.score_raw = -9999;

cdata.GT = [];
%% Find initial matches 
[ matchInfo ]= make_initialmatches_mcho2( cdata.view, cdata.mparam, 'verbose', true );
%[ matchInfo ]= make_initialmatches_mcho( cdata.view, cdata.mparam, 'verbose', true );
if 0
    hFigIM = figure('NumberTitle','off'); 
    set(hFigIM, 'Name',['Initial Matching: ' num2str(numel(matchInfo)) ]);
    cdata.bPair = 1;    cdata.matchInfo = matchInfo;
    cdata.nInitialMatches = size(matchInfo.match,2);
    showInitialMatches2;
    pause;
end
drawnow;
cand_matchlist = matchInfo.match;
cand_matchdist = matchInfo.dist;

[ uniq_feat1, tmp, new_feat1 ] = unique(cand_matchlist(1,:));    
[ uniq_feat2, tmp, new_feat2 ] = unique(cand_matchlist(2,:));
matchlist_tmp = [ new_feat1; new_feat2 ];
%matchlist_tmp = cand_matchdist;
if 1
    %fprintf('- initial candidates: %d \n\n',size(cand_matchlist,1));
    edgeAttr1 = computeEdgeAttr( cdata.view(1).frame(:,uniq_feat1) );
    tic
    edgeAttr2 = computeEdgeAttr( cdata.view(2).frame(:,uniq_feat2) );
    toc
    featBin1 = makeFeatBin(edgeAttr1);
    tic
    featBin2 = makeFeatBin(edgeAttr2);
    toc
    tic
    [ eSimVal ] = computeDotProdSimilarity_sym( int32(matchlist_tmp), featBin1, featBin2); % symmetric affinities
    %HoughInitial;
    %wE = ones(size(edgeAttr1));
    %wE(:,:,1) = 0.1;    wE(:,:,2) = 0.1;    wE(:,:,3) = 0.05;
    %[ eSimMat0 eSimVal ] = computeSimilarity( int32(matchlist_tmp), wE, edgeAttr1, edgeAttr2 );
    %cdata.affinityMatrix = computeSimilarity2( int32(matchlist_tmp), wE, eSimVal); 
    cdata.affinityMatrix = double(eSimVal);
    toc
else
    eSimVal = [];
    [ distance tmp] = affineTransferDistanceMEX_R( int32(matchlist_tmp'), cdata.view(1).frame(:,uniq_feat1), cdata.view(2).frame(:,uniq_feat2), 0, 0);
    %[ distance tmp] = affineTransferDistanceMEX_R( int32(cand_matchlist), cdata.view(1).affMatrix(:,:), cdata.view(2).affMatrix(:,:), 0, 0);
    cdata.affinityMatrix = max(100 - distance, 0);
end
%% Make the overlapping groups of initial matches
%tic
[ cdata.group1 cdata.group2 ] = make_group12(matchlist_tmp(1:2,:));
[ cdata.L12, cdata.E12, cdata.nP1, cdata.nP2 ] = updateL12E12( matchlist_tmp' );
%toc

unaryWeight = 10.0*ones(size(cand_matchdist));
vSimVal = max( 0, 0.8 - cand_matchdist );
try
    cdata.affinityMatrix(1:(size(cdata.affinityMatrix,1)+1):end)...
        = unaryWeight(matchlist_tmp(1,:)) .* vSimVal;
catch
    pause;
end

if 1%aparam.bOneToOne 
    % eliminate conflicting elements to prevent conflicting walks
    cdata.affinityMatrix = cdata.affinityMatrix.*~getConflictMatrix(cdata.group1, cdata.group2);
end
%cdata.affinityMatrix(1:(size(cdata.affinityMatrix,1)+1):end) = 0;

%% select elements of the tree structure
% lookup table 
% nNode_model = size(cdata.view(1).feat,1);
% star
if 0
    validMatrix = zeros(size(cdata.affinityMatrix));
    for i = 1:size(model.graph_m.edges_star,1)
        v1 = find( cand_matchlist(:,1) == model.graph_m.edges_star(i,1));
        v2 = find( cand_matchlist(:,1) == model.graph_m.edges_star(i,2));
        validMatrix( v1, v2) = true;
        validMatrix( v2, v1) = true;
    end
    cdata.affinityMatrix = cdata.affinityMatrix.*validMatrix;
end
fprintf('matrix density: valid / full = %d / %d (%5.2f)\n', nnz(cdata.affinityMatrix), prod(size(cdata.affinityMatrix)),...
    nnz(cdata.affinityMatrix) / prod(size(cdata.affinityMatrix)));

%% ------------------ Perform graph matching
tic;
[ X_raw ] = eval(strFunc);
res.time = toc;
if method.postDiscretization 
    X_sol = greedyMapping(X_raw, cdata.group1, cdata.group2);
    score_raw = X_sol'*cdata.affinityMatrix*X_sol;
else
    X_sol = X_raw;
end
%sum(X_sol - X_raw)
score_GM = X_sol'*cdata.affinityMatrix*X_sol;
%score_GM
%fprintf('score_raw: %f score_GM: %f\n', score_raw, score_GM);


%% ----------------- Evaluate the solutions
%X_GT = extrapolateGT( cdata.view, cand_matchlist, cdata.GT,  cdata.mparam.extrapolation_dist ); % extrapolate the groundtruths
%X_sol_EXT = extrapolateMatchIndicator( cdata.view, cand_matchlist, X_sol, cdata.mparam.extrapolation_dist ); % extrapolate the solutions

matchIdx_GM = find(X_sol);
matchScore_GM = X_raw(matchIdx_GM);     matchScore_GM = matchScore_GM./sum(matchScore_GM);
%matchList_GM = cand_matchlist(matchIdx_GM,:);
%nDetected = nnz(X_sol_EXT);    nTrue = nnz(X_GT);  nTP = nnz(X_GT & X_sol_EXT );
%recall_GM = nTP/nTrue;
%perform_data(iterGM,1:7) = [ size(cand_matchlist,1), nTrue, nDetected, nTP, score_raw, score_GM, res.time];
%if iterGM == 1
%    score_GM_oneshot = score_GM;        nTP_GM_oneshot = nTP;
%end
%inlierness_GM = sum(cdata.affinityMatrix(find(X_sol),find(X_sol)),1);
%cand_matchlist(find(X_sol),:)'
%sum(cdata.affinityMatrix(find(X_sol),find(X_sol)),1)

    res.X = X_sol;
    res.X_raw = X_raw;
    res.score_raw = score_raw;
    res.score_GM = score_GM;
    res.cand_matchlist = cand_matchlist;
    res.affinityMatrix = cdata.affinityMatrix;
    res.eSimVal = eSimVal;
    res.vSimVal = vSimVal;
%end

%fprintf('result of iter #%d returned\n', iterGM);
%res.perform_data = perform_data(1:min(iterGM+1,maxIterGM),:);
