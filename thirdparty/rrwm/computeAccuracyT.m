function pY = computeAccuracyT( dsave, distRatioTrue )

pY = zeros(1,100); 
nImage = 0;
for i=1:numel(dsave)
    
    %nGT
    nGT = size(dsave(i).ptGT,2);
    if nGT < 2
        continue;
    end
    
    nImage = nImage + 1;        
    [ tmp, idbest ] = sort(dsave(i).X_raw,'descend');
    
    %trueGT = zeros(1,size(dsave(i).ptGT,2));
    nMatch = size(dsave(i).match,2);
    ptGT = dsave(i).ptGT(:,dsave(i).match(1,:));
    ptM = dsave(i).ptQ(:,dsave(i).match(2,:));
    endpointErr = sqrt(sum((ptGT - ptM).^2,1))/sqrt(sum(dsave(i).object_size.^2));
    bTrue = endpointErr < distRatioTrue;
    
    sorted_true = bTrue(idbest);
    sorted_id = dsave(i).match(1,idbest);
%     if numel(sorted_true) < maxMn
%         sorted_true(maxMn) = 0;
%     else
%         sorted_true = sorted_true(1:maxMn);
%     end
    
    for j=1:100
        nJ = ceil(j*nMatch/1000);
        pY(j) = pY(j) + numel(unique(sorted_id(sorted_true(1:nJ))))/nGT;
    end
    %pY = cumsum(sorted_true);
    
    
end

pY = pY / nImage;

end