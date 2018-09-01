function [ exmp ] = makeExampleSet_img( classid, dataset )
% make example sets for traing a predictor
idxValid = find([ dataset.classid ] == classid );
idxValid = idxValid(randperm(numel(idxValid))); % random permutation!

for i=1:numel(idxValid)
    exmp(i).type = dataset(idxValid(i)).type;
    exmp(i).classid = dataset(idxValid(i)).classid;         % class label
    exmp(i).dataid =  dataset(idxValid(i)).dataid;
    view =  extract_localfeatures_wrapper( dataset(idxValid(i)).imgpath, dataset(idxValid(i)).tmppath );
    exmp(i).ptGT = dataset(idxValid(i)).ptGT;
    idxSel = [];
    pt = exmp(i).ptGT;
    for j=1:size(pt,2)
        idx1 = nearestneighbour(pt(:,j), view.frame(1:2,:), 'NumberOfNeighbours', 1);
        idxSel = [idxSel idx1];
    end
    exmp(i).seqGT = idxSel;
    exmp(i).frameGT = view.frame(:,idxSel);
    
    exmp(i).fname = dataset(idxValid(i)).imgpath;
    exmp(i).tname = dataset(idxValid(i)).tmppath;
    exmp(i).box = [];%dataset(idxValid(i)).bbox;
    exmp(i).size = [ size(view.img,1) size(view.img,2) ];
    exmp(i).object_size = [ max(dataset(idxValid(i)).ptGT(1,:)) max(dataset(idxValid(i)).ptGT(2,:)) ]...
        - [ min(dataset(idxValid(i)).ptGT(1,:)) min(dataset(idxValid(i)).ptGT(2,:)) ];
    %exmp(i).label = 2*( classid == dataset(idxValid(i)).classid ) - 1;
    %exmp(i).fV = [];
    %exmp(i).fE = [];
    %exmp(i).score = 0;
end