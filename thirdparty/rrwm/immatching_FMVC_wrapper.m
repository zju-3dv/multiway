function [ res cdata ] = immatching_FMVC_wrapper( x, view_m, svm, y )
bVerbose = false;

iter = 1;
idxAlg = 1;
set_alg;
set_param_GM;

cdata.fparam = fparam;  cdata.aparam = aparam;  cdata.mparam = mparam;

view_m.desc = svm.wV;
fE1 = svm.wE;
cdata.view(1) = view_m;

% command to call a matching function
strFunc = ['feval(@' func2str(alg(idxAlg).fhandle)];
for j = 1:length(alg(1).variable)
    strFunc = [ strFunc ', cdata.' alg(idxAlg).variable{j} ];
end
strFunc = [strFunc ')'];

cdata.view(2) =  extract_localfeatures_wrapper( x.fname, x.tname );
        
[ matchInfo ]= make_initialmatches_mcho2( cdata.view, cdata.mparam );
if 0%bShow
    hFigIM = figure('NumberTitle','off');   iptsetpref('ImshowBorder','tight');
    %set(hFigIM, 'Name',['Initial Matching: ' num2str(size(matchInfo.match,2) ]);
    showInitialMatches2;
    hold on;
end

%fprintf('- initial candidates: %d \n',size(matchInfo.match,2));
[ uniq_feat2, tmp, new_feat2 ] = unique(matchInfo.match(2,:));
cand_matchlist_uniq = [ matchInfo.match(1,:); new_feat2' ];        
edgeAttr2 = computeEdgeAttr( cdata.view(2).frame(:,uniq_feat2) );
%fE2 = svm.maskE(1) * makeFeatBin(edgeAttr2);
fE2 = makeFeatBin(edgeAttr2);
[ eSimVal ] = computeDotProdSimilarity_sym( int32(cand_matchlist_uniq), fE1, fE2); % symmetric affinities


%% margin rescaling 
pt1 = y.pt(:,matchInfo.match(1,:));
pt2 = cdata.view(2).frame(1:2,matchInfo.match(2,:));
lossScore = sum((pt1 - pt2).^2,1) > y.distThres^2;
%lossScore = sum(min(1.0, sum((pt1 - pt2).^2,1)./y.distThres^2),1);
lossScore = lossScore / size(svm.wV,2);
%%

if 0
    vSimVal = svm.maskV * matchInfo.sim + lossScore;
    % make the affinity matrix for GM
    cdata.affinityMatrix = double(eSimVal);
    cdata.affinityMatrix(1:(size(cdata.affinityMatrix,1)+1):end) = vSimVal;
    %pause;
else
    %vSimVal = svm.maskV * matchInfo.sim + 1.0*lossScore;
    vSimVal = matchInfo.sim + lossScore;
    mat_vSimVal = repmat(vSimVal, [ numel(vSimVal) 1 ])/(2*(size(svm.wV,2)-1)); % distribute the v score
    eSimVal = eSimVal + mat_vSimVal + mat_vSimVal';
    eSimVal(1:size(eSimVal,1)+1:end) = 0;
    minSimE = min(eSimVal(:));
    if minSimE<0
        %eSimVal(eSimVal<0) = 0;
        eSimVal = eSimVal - minSimE; % make them non-negative
        %eSimVal(eSimVal~=0) = eSimVal(eSimVal~=0) - minSimE; % make them non-negative
    end
    % make the affinity matrix for GM
    cdata.affinityMatrix = double(eSimVal);
    cdata.affinityMatrix(1:(size(cdata.affinityMatrix,1)+1):end) = 0;
    %pause; 
end

%% check the overlapping groups of initial matches
[ cdata.group1 cdata.group2 ] = make_group12(matchInfo.match(1:2,:));
if 1%aparam.bOneToOne 
    % eliminate conflicting elements to prevent conflicting walks
    cdata.affinityMatrix = cdata.affinityMatrix.*~getConflictMatrix(cdata.group1, cdata.group2);
end

if bVerbose 
    fprintf('matrix density: valid / full = %d / %d (%5.2f)\n', nnz(cdata.affinityMatrix), prod(size(cdata.affinityMatrix)),...
    nnz(cdata.affinityMatrix) / prod(size(cdata.affinityMatrix)));
end
%% ------------------ Perform graph matching
tic;
[ X_raw ] = eval(strFunc);
res.time = toc;
if alg(idxAlg).postDiscretization 
    X_sol = greedyMapping(X_raw, cdata.group1, cdata.group2);
else
    X_sol = X_raw;
end

%sum(X_sol - X_raw)
score_GM = X_sol'*cdata.affinityMatrix*X_sol;
if bVerbose, fprintf('* iter #%d - score_GM: %f\n', iter, score_GM); end
matchIdx_GM = find(X_sol);
matchScore_GM = X_raw(matchIdx_GM);     matchScore_GM = matchScore_GM./sum(matchScore_GM);

%idxframe = matchInfo.match(2,find(X_sol));    % bug here
idxframe(matchInfo.match(1,find(X_sol))) = matchInfo.match(2,find(X_sol));    % be careful!!!!!
fV = cdata.view(2).desc(:,idxframe);
edgeAttr_t = computeEdgeAttr( cdata.view(2).frame(:,idxframe) );
fE = makeFeatBin(edgeAttr_t);  
pt = cdata.view(2).frame(1:2,idxframe); 
sV = sum(sum(svm.wV.*fV));
sE = sum(sum(svm.wE.*fE));


res.X = X_sol;
res.X_raw = X_raw;
res.score_GM = score_GM;
res.fE = fE;
res.fV = fV;
res.sV = sV;
res.sE = sE;
res.pt = pt;
res.seq = idxframe;
res.frame = cdata.view(2).frame;
%res.desc = cdata.view(2).desc;
res.cand_matchlist = matchInfo.match;
res.affinityMatrix = cdata.affinityMatrix;
end