function [ E ] = makeExampleSet( classid, dataset )
% make example sets for traing a binary classifier

for i=1:numel(dataset)
    E(i).type = dataset(i).type;
    E(i).fname = dataset(i).imgpath;
    E(i).tname = dataset(i).tmppath;
    E(i).box = [];%dataset(i).bbox;
    E(i).label = 2*( classid == dataset(i).classid ) - 1;
    E(i).fV = [];
    E(i).fE = [];
    E(i).score = 0;
end