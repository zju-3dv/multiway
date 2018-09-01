function matchInfo = initMatch(viewPair,kNN,nMax,thScore,thRatio,thDist)

% Find k nearest-neighbor matches between feature sets
% matchInfo(:).match = [ featA featB ]
% matchInfo(:).sim  = similarity between featA & featB

% global idx
% desc1 = normc(viewPair(1).desc(idx,:));
% desc2 = normc(viewPair(2).desc(idx,:));
desc1 = normc(viewPair(1).desc);
desc2 = normc(viewPair(2).desc);

match = [];
sim = [];

% n = min(kNN,size(desc2,2));
% for i = 1:size(desc1,2)
%     dotprods = single(desc1(:,i))'*single(desc2);
%     [vals,indx] = sort(dotprods,'descend');
%     if vals(1)/(vals(2)+eps) >= thRatio
%         k = sum( vals(1:n) >= thScore );
%         match = [match [i*ones(1,k);indx(1:k)]];
%         sim = [sim vals(1:k)];
%     end
% end

frame2 = viewPair(2).frame(1:2,:);
thDist = max(size(viewPair(2).img))*thDist;
for i = 1:size(desc1,2)
    dotprods = single(desc1(:,i))'*single(desc2);
    [vals,indx] = sort(dotprods,'descend');
    if vals(1)/(vals(2)+eps) >= thRatio
        c = 0;
        for j =1:length(indx)
            if vals(j) > thScore
                c = c+1;
                match = [match [i;indx(j)]];
                sim = [sim vals(j)];
                dist = sqrt(sum(bsxfun(@minus,frame2(:,indx),frame2(:,indx(j))).^2,1));
                vals(dist<thDist) = -inf;
                if c >= kNN
                    break
                end
            end
        end
    end
end

if length(sim) > nMax
    [vals,indx] = sort(sim,'descend');
    sim = vals(1:nMax);
    match = match(:,indx(1:nMax));
end

matchInfo.sim = single(sim);
matchInfo.match = uint16(match);
    


