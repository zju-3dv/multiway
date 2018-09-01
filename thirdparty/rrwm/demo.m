clear
addpath gm-toolbox

load data
img1 = obj1.img; % images
img2 = obj2.img;
desc1 = obj1.desc; % descriptors
desc2 = obj2.desc;
xy1 = obj1.frame(1:2,:); % feature locations
xy2 = obj2.frame(1:2,:);

wEdge = 10; % weigth of rigidicity
% parameter for initial candidate search
kNNInit = 20; % number of candidate correpondences for each feature
nMaxInit = 1000; % max number of candidate correpondences in total
thScore = 0.5; % affinity threshold below which a pair will not be considered
thRatio = 1; % 

%% obtain initial matches

match = [];
sim = [];
n = min(kNNInit,size(desc2,2));

for i = 1:size(desc1,2)
    dotprods = single(desc1(:,i))'*single(desc2);
    [vals,indx] = sort(dotprods,'descend');
    if vals(1)/(vals(2)+eps) >= thRatio
        k = sum( vals(1:n) >= thScore );
        match = [match [i*ones(1,k);indx(1:k)]];
        sim = [sim vals(1:k)];
    end
end

if length(sim) > nMaxInit
    [vals,indx] = sort(sim,'descend');
    sim = vals(1:nMaxInit);
    match = match(:,indx(1:nMaxInit));
end

%% solve graph matching
if wEdge <= 0
    Xraw = double(sim);
else
    % construct graphs
    fprintf('- constructing graphs ...\n');
    [uniq_feat1,~,new_feat1] = unique(match(1,:));
    [uniq_feat2,~,new_feat2] = unique(match(2,:));
    cand_matchlist_uniq = [new_feat1';new_feat2'];
    edgeAttr1 = computeEdgeAttr(xy1(:,uniq_feat1));
    edgeAttr2 = computeEdgeAttr(xy2(:,uniq_feat2));
    fE1 = makeFeatBin(edgeAttr1);
    fE2 = makeFeatBin(edgeAttr2);
    eSimVal = computeDotProdSimilarity_sym( int32(cand_matchlist_uniq), fE1, fE2); % symmetric affinities
    vSimVal = sim;
    mSimVal = repmat(vSimVal,[numel(vSimVal),1]); % distribute the v score
    affinityMatrix = wEdge*eSimVal + mSimVal + mSimVal';
    affinityMatrix(1:size(affinityMatrix,1)+1:end) = 0;
    affinityMatrix(affinityMatrix<1e-3) = 0;
    affinityMatrix = sparse(double(affinityMatrix));
    fprintf('- running random walk matching ...\n');
    [group1,group2] = make_group12(match(1:2,:));
    Xraw = RRWM(affinityMatrix,group1,group2);
end

%% discretization
X = greedyMatch(match,Xraw);

figure;
img(1:size(img1,1),1:size(img1,2),:) = img1;
img(1:size(img2,1),1+size(img1,2):size(img2,2)+size(img1,2),:) = img2;
imshow(img); hold on
xy2(1,:) = xy2(1,:) + size(img1,2);
for i = 1:size(match,2)
    if X(i)
        plot([ xy1(1,match(1,i)), xy2(1,match(2,i)) ]...
            ,[ xy1(2,match(1,i)), xy2(2,match(2,i)) ],...
            '-+','LineWidth',1,'MarkerSize',5,...
            'color', 'y');
    end
end



