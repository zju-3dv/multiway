function [ voting nAddedVote ] = voteCandidate( all_feat1, all_type1, all_feat2, all_type2, matchList_GM, scoreAnchor, anchorFeat1, Ti21, k_neighbor1, k_neighbor2  )

bDebugMode = 0;
lambda_anchor = 0.0; 
lambda_dist = 1.0;
% bounded search from all feat1 of the same type
ptAnchor1 = all_feat1(anchorFeat1,1:2);
typeAnchor1 = all_type1(anchorFeat1);
  
ind_cand0 = find( all_type1 == typeAnchor1);
if length(ind_cand0) > k_neighbor1%+1
    %ind_featCand1_tmp=BruteSearchMex(all_feat1(ind_cand0,1:2)',ptAnchor1','r',radius_neighbor1);
    ind_featCand1_tmp=BruteSearchMex(all_feat1(ind_cand0,1:2)',ptAnchor1','k',k_neighbor1);%+1);
else
    ind_featCand1_tmp=BruteSearchMex(all_feat1(ind_cand0,1:2)',ptAnchor1','k',length(ind_cand0));
end
ind_featCand1 = ind_cand0(ind_featCand1_tmp);
%ind_featCand1 = setdiff(ind_featCand1,anchorFeat1); % ignore itself

if bDebugMode
    plot(all_feat1(ind_featCand1,1),all_feat1(ind_featCand1,2),'o','LineWidth',1,'MarkerSize',7,'MarkerFaceColor','m','color','w');
end
% project the neighbors onto img2
tmp1 = Ti21 * [ all_feat1(ind_featCand1,1:2) ones(length(ind_featCand1),1) ]';
ptAnchor1_projected = tmp1(1:2,:)';
if bDebugMode
    plot(ptAnchor1_projected(:,1)+img_offset,ptAnchor1_projected(:,2),'x','LineWidth',2,'MarkerSize',7,'MarkerFaceColor','m','color','w');
end

% make a lookup table of feat2 matching feat1 in matches of GM
ind_feat2_matching_feat1 = zeros(size(all_feat1,1),1);
ind_feat2_matching_feat1( matchList_GM(:,1), 1) = matchList_GM(:,2);

voting = sparse(size(all_feat1,1),size(all_feat2,1));
%voting = zeros(size(all_feat1,1),size(all_feat2,1));
nAddedVote = 0;
for iter_j = 1:length(ind_featCand1)
    ptCand1 = all_feat1(ind_featCand1(iter_j),1:2);
    typeCand1 = all_type1(ind_featCand1(iter_j));
    ptCand1_projected = ptAnchor1_projected(iter_j,:);
    % draw the match line
    if bDebugMode
        plot([ ptCand1(1), ptCand1_projected(1)+img_offset ],[ ptCand1(2), ptCand1_projected(2) ],...
            '-','LineWidth',2,'MarkerSize',10,'color', 'b');
        circle_xyr(ptCand1_projected(1)+img_offset,ptCand1_projected(2),radius_neighbor2);
    end
    % bounded search from all feat2 of the same type
    %ind_cand_feat2 = find( all_type2 > -1);
    ind_cand_feat2 = find( all_type2 == typeCand1);
    if length(ind_cand_feat2) > k_neighbor2
        %[ ind_featCand2 dist_featCand2 ] = BruteSearchMex(all_feat2(:,1:2)',ptCand1_projected','r',radius_neighbor2);
        [ ind_featCand2_tmp dist_featCand2 ] = BruteSearchMex(all_feat2(ind_cand_feat2,1:2)',ptCand1_projected','k',k_neighbor2);
    else
        [ ind_featCand2_tmp dist_featCand2 ] = BruteSearchMex(all_feat2(ind_cand_feat2,1:2)',ptCand1_projected','k',size(all_feat2,1));
    end
    ind_featCand2 = ind_cand_feat2(ind_featCand2_tmp);
        
    if bDebugMode
        featCand2 = all_feat2(ind_featCand2,:); 
        plot(featCand2(:,1)+img_offset,featCand2(:,2),'o','LineWidth',1,'MarkerSize',7,'MarkerFaceColor','b','color','w');
    end
    %dist_featCand2
    likelihood_cand = exp(-lambda_dist*dist_featCand2); 
    %pause;
    iind_match_of_GM = find(ind_featCand2 == ind_feat2_matching_feat1(ind_featCand1(iter_j)));
    [ tmp, iind_max ] = max(likelihood_cand);
    % if the voting includes the match of the current GM
    if ~isempty(iind_match_of_GM) && ( iind_match_of_GM == iind_max )
        ind_featCand2 = ind_featCand2(iind_match_of_GM);
        likelihood_cand = 1; % obtain all the votes
    end
    if isempty(ind_featCand2), continue; end
    likelihood_cand=likelihood_cand./sum(likelihood_cand);
    sim_cand = scoreAnchor^lambda_anchor*likelihood_cand;
    voting(ind_featCand1(iter_j), ind_featCand2) = voting(ind_featCand1(iter_j), ind_featCand2) + sim_cand;
    %voting(ind_featCand1(iter_j), ind_featCand2(iter_k)) = max( voting(ind_featCand1(iter_j), ind_featCand2(iter_k)), sim_cand);
    nAddedVote = nAddedVote + length(ind_featCand2);
    if bDebugMode
        %drawEllipse3( featCand(iter_k,1:5) + [ img_offset 0 0 0 0], 1,'b',1);
        pause;
    end
end