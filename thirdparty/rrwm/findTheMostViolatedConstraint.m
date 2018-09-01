function [ res ] = findTheMostViolatedConstraint( w, x, y )
% find the most violated constraint in the current GM
% based on
% margin rescaling: argmax_y delta(yi, y) + <psi(x,y), w>
% slack rescaling: argmax_y delta(yi, y) (1 + <psi(x,y), w> - <psi(x,yi), w>)

set_alg;

nP1 = numel(y.seq); % # of pts in model
nP2 = size(x.view.frame,2);  % # of pts in test input

% extract features from input
%tic
edgeAttr_t = computeEdgeAttr( x.view.frame );
featBin_t = makeFeatBin(edgeAttr_t);
%toc

% reshape the weight w
wE = single(reshape(w, size(featBin_t,1), []));

E12 = ones(nP1,nP2);
[L12(:,1) L12(:,2)] = find(E12);     % match list
[group1 group2] = make_group12(L12'); % group info
%fprintf('-# of initial candidates: %d \n',size(cand_matchlist,1));

% compute the affinity matrix for GM
[ eSimVal ] = computeDotProdSimilarity_sym( int32(L12'), wE, featBin_t); % symmetric affinities
%eSimVal = eSimVal - min(eSimVal(:)); % make them non-negative

%delta = (y_seq == ybar_seq)/numel(y_seq); 
loss_score = zeros(nP1, nP2);
for i = 1:nP1, loss_score(i,y.seq(i)) = 1; end
%loss_score = loss_score'; % wrong
loss_score = (1 - loss_score(:))/nP1;


% make affinity matrix
if 1
    % margin rescaling: argmax_y delta(yi, y) + (<psi(x,y), w> - <psi(x,yi), w>)
    vSimVal = 0.5*loss_score;
    %M = double(eSimVal);
    %M(1:(size(M,1)+1):end) = vSimVal;
    mat_vSimVal = repmat(vSimVal, [ 1 numel(vSimVal) ])/(2*(nP1-1)); % distribute the v score
    eSimVal = eSimVal + mat_vSimVal + mat_vSimVal';
    eSimVal(1:size(eSimVal,1)+1:end) = 0;
    M = double(eSimVal);
else
    % slack rescaling: argmax_y delta(yi, y) (1 + <psi(x,y), w> - <psi(x,yi), w>)
    % ??
end
%size(M)

%% Ground Truth
GT.seq = y.seq; GT.matrix = zeros(nP1, nP2);
for ip = 1:nP1, GT.matrix(ip,y.seq(ip)) = true; end
GT.bool = GT.matrix(:);

%% params for GM
problem.nP1 = nP1;  problem.nP2 = nP2;
%problem.P1 = model.view_m.frame(1:2,:)';
%problem.P2 = x.view.frame(1:2,:)';
problem.L12 = L12;  problem.E12 = E12;
problem.group1 = group1;    problem.group2 = group2;
problem.GTbool = GT.bool;   problem.GTseq = GT.seq;
problem.affinityMatrix = M;

for j = 1:1%length(method_GM)
    [accuracy matchScore time X] = wrapper_GM(alg(j), problem);
end

if nnz(X) == numel(y.seq)
    [ ind1 ind2 ] = find(reshape(X,nP1,nP2)'); %correct
    %[ ind1 ind2 ] = find(reshape(X,nP2,nP1)); %wrong
    res.seq = ind1';
    res.pt = x.view.frame(1:2,res.seq);
else
    rand_y_seg = randperm(nP2);
    res.seq = rand_y_seg(1:numel(y.seq));
    res.pt = x.view.frame(1:2,yhat.seq);
end
    
    