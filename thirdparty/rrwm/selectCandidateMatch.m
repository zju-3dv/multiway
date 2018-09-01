function [ cand_matchlist cand_matchdist ] = selectCandidateMatch( voting_space, feat1, desc1, feat2, desc2, max_candidates, threshold_dissim, mparam )
% collect new candidate matches from the voting space
% : making sure that each feature has at most k matches
if nargin < 8
    mparam.kNN = 10;
    mparam.bFilterMatch = 1;                   % filter out overlapping matches (very similar, redundant matchs)
    mparam.redundancy_thres = 2.0;             % criterion to filter out overlapping matches (pixel distance)
end
[ ind_feat1, ind_feat2, vote_score] = find(voting_space);
[ tmp, iv ] = sort(vote_score,'descend');
nSelMatch = 0;   iter_j = 0; sel_iv=[]; sel_sqdist = [];
nMax = length(iv);
if nMax < 1
    cand_matchlist = [];
    return;
end
mparam.bFilterMatch = 0;
numMatch1 = zeros(size(voting_space,1),1);
numMatch2 = zeros(size(voting_space,2),1);

if mparam.bFilterMatch
    bSelf = 0;
    %tmp = find(isinf(tmp),1,'last');
    %if isempty(tmp)
        nSelMatch = 0;
    %else
    %    nSelMatch = min(tmp,nMax);
    %end
    %sel_iv = iv(1:nSelMatch); % find first matches with Inf values
    XY_sel_view1 = [];%(1:nSelMatch,:) = feat1(ind_feat1(sel_iv),1:2);
    XY_sel_view2 = [];%(1:nSelMatch,:) = feat2(ind_feat2(sel_iv),1:2);
    for  iter_cand = 1:nMax%nSelMatch+1:nMax
        cur_iv = iv(iter_cand);
        if nSelMatch >= max_candidates, break; end
        XY_view1 = feat1(ind_feat1(cur_iv),1:2);
        XY_view2 = feat2(ind_feat2(cur_iv),1:2);
        
        if iter_cand > 1
            if ~bSelf
                %finds the points within the search radius
                ridx_view1=BruteSearchMex(XY_sel_view1',XY_view1','r',mparam.redundancy_thres);
                ridx_view2=BruteSearchMex(XY_sel_view2',XY_view2','r',mparam.redundancy_thres);
            else
                %finds the points within the search radius
                ridx_view1=BruteSearchMex([XY_sel_view1; XY_sel_view2]',XY_view1','r',mparam.redundancy_thres);
                ridx_view2=BruteSearchMex([XY_sel_view2; XY_sel_view1]',XY_view2','r',mparam.redundancy_thres);
            end
            equi_ridx = intersect(ridx_view1, ridx_view2);
        else
            equi_ridx = [];
        end
        
        if isempty(equi_ridx) % if there's no equi matches
            desc_sqdist = sum((desc1(ind_feat1(cur_iv),:)-desc2(ind_feat2(cur_iv),:)).^2);
            if desc_sqdist <= threshold_dissim^2
                nSelMatch = nSelMatch + 1;
                % insert it into the selected xy list
                sel_iv = [ sel_iv; cur_iv];
                sel_sqdist = [sel_sqdist; desc_sqdist ]; 
                XY_sel_view1(nSelMatch,:) = XY_view1;
                XY_sel_view2(nSelMatch,:) = XY_view2;
            end
        end            
    end
else
    while nSelMatch < max_candidates && iter_j < length(iv)
        iter_j = iter_j + 1;
        cur_iv = iv(iter_j);
        if numMatch1(ind_feat1(cur_iv)) < mparam.kNN %&& numMatch2(ind_feat2(cur_iv)) < mparam.kNN
            %if cdata.fparam.bFeatExtUse(cdata.view(1).typeFeat(ind_feat1(cur_iv)))
                desc_sqdist = sum((desc1(ind_feat1(cur_iv),:)-desc2(ind_feat2(cur_iv),:)).^2);
                if desc_sqdist <= threshold_dissim^2
                    numMatch1(ind_feat1(cur_iv)) = numMatch1(ind_feat1(cur_iv)) + 1;
                    numMatch2(ind_feat2(cur_iv)) = numMatch2(ind_feat2(cur_iv)) + 1;
                    nSelMatch = nSelMatch + 1;
                    sel_iv = [ sel_iv; cur_iv];
                    sel_sqdist = [sel_sqdist; desc_sqdist ];                 
                end
            %end
        end
    end
end
%sel_iv = iv(1:min(max_candidates,length(iv)));
cand_matchlist = [ ind_feat1(sel_iv) ind_feat2(sel_iv) ]; 
cand_matchdist = sel_sqdist;