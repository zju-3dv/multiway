function pMatch = runGraphMatchBatch(datapath,viewList,pairingMode,pairingParam,varargin)

% Options:
% pairingMode = 'all' | 'loop' | 'neighbor' | 'userinput'
% pairingParam = ...
% the radius of neighnorhood, if pairingMode == 'loop' and 'neighbor'
% a binary matrix indicating the pairs to match, if pairingMode == 'userinput'

nImg = length(viewList);
pairs = false(nImg,nImg);
switch lower(pairingMode)
    case 'all'
        pairs = true(nImg,nImg);
    case 'loop'
        for i = 1:size(pairs,1)
            idx = i+1:i+pairingParam;
            idx(idx>nImg) = mod(idx(idx>nImg),nImg);
            pairs(i,idx) = true;
        end
        pairs = pairs | pairs';
    case 'neighbor'
        for i = 1:size(pairs,1)
            idx = i+1:i+pairingParam;
            idx(idx>nImg) = [];
            pairs(i,idx) = true;
        end
        pairs = pairs | pairs';
    case 'userinput'
        pairs = pairingParam;
    otherwise
        error('pairingMode invalid!!');
end

for i = 1:length(viewList)
    views(i) = load(sprintf('%s/%s',datapath,viewList{i}));
end

pMatch(length(viewList),length(viewList)) = struct('matchInfo',[],'nFeature',[],'Xraw',[],'X',[],'filename',[]);

for i = 1:length(viewList)
    for j = i+1:length(viewList)
        if pairs(i,j)
            fprintf('- running graph matching ...\n');
            mData = graphMatch(views([i,j]),varargin{:});
            mData.filename = viewList([i,j]);
            fprintf('- Done! (%d,%d) feautres, %d initial matches, %d output matches\n',...
                mData.nFeature(1),mData.nFeature(2),length(mData.X),sum(mData.X));
            pMatch(i,j) = mData;
        end
    end
end
end