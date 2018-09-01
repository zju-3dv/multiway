function [ idx1 ID1 idx2 ID2 dumVal dumSize] = make_groups_slack( idx1, ID1, idx2, ID2 )

nG1 = ID1(end);
nG2 = ID2(end);
dumVal = nG2 - nG1;
dumSize = nG2;

maxIdx = length(idx1);
addIdx = ((maxIdx+1):(maxIdx+nG2))';

idx1 = [idx1; addIdx];
ID1 = [ID1; (nG1+1)*ones(nG2,1)];

addID = (1:nG2)';
tempIdx = zeros(maxIdx+nG2,1);
tempID = zeros(maxIdx+nG2,1);

i = 1;
j = 0;
while i < maxIdx
    tempIdx(i+j) = idx2(i);
    tempID(i+j) = ID2(i);
    if ID2(i) ~= ID2(i+1)
        j = j+1;
        tempIdx(i+j) = addIdx(j);
        tempID(i+j) = addID(j);
    end
    i = i+1;
end

while i <= maxIdx
    tempIdx(i+j) = idx2(i);
    tempID(i+j) = ID2(i);
    if i == maxIdx
        j = j+1;
        tempIdx(i+j) = addIdx(j);
        tempID(i+j) = addID(j);
    end
    i = i+1;
end

idx2 = tempIdx;
ID2 = tempID;

end


