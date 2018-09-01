function [ matchInfo ] = make_initialmatches_mcho2( viewInfo, mparam, varargin )
% Find NN matches between feature sets
% matchInfo(:).match = [ featA featB ]
% matchInfo(:).dist  = descriptor distance between featA & featB;
%
% Minsu Cho, Seoul National University
% Updated: 16th November March 2011

bVerbose = false;
bShow = false;

for k=1:2:length(varargin)
  opt=lower(varargin{k}) ;
  arg=varargin{k+1} ;
  switch opt
    case 'verbose'
      bVerbose = arg ;
    case 'showmatches'
      bShow = arg ;
    otherwise
      error(sprintf('Unknown option ''%s''', opt)) ;
  end
end

initialMatch = [];
simInitialMatch = [];

if length(viewInfo) == 1
    bSelf = 1;    view1 = 1;  view2 = 1;
else
    bSelf = 0;    view1 = 1;  view2 = 2;
end

% set feature types to match 
if any(strcmp(fieldnames(mparam), 'bFeatExtUse')) 
    typeToMatch =  mparam.bFeatExtUse;
else
    typeToMatch = [ viewInfo(view1).nFeatOfExt > 0 ] & [ viewInfo(view2).nFeatOfExt > 0];
end

maxFeatType = max(viewInfo(view1).type);
nMatchOfFeatType = zeros(maxFeatType,1);
nFeatTypeToMatch = nnz(typeToMatch);%length(unique(viewInfo(1).typeFeat));

for i=1:maxFeatType
    
    if ~typeToMatch(i),  continue;   end

    range1 = find( viewInfo(view1).type == i );
    range2 = find( viewInfo(view2).type == i );    
    
    if isempty(range1) || isempty(range2),  continue;   end
    
    tic;
    if mparam.bReflective
        [ tmpInitialMatch, simdot ] = descmatch_dot( viewInfo(view1).desc(:,range1),...
            viewInfo(view2).desc_ref(:,range2),bSelf , mparam.distRatio, mparam.distThres, mparam.kNN, 1 ) ;
    else
        %try
        [ tmpInitialMatch, simdot ] = descmatch_dot( viewInfo(view1).desc(:,range1),...
            viewInfo(view2).desc(:,range2), bSelf , mparam.distRatio, mparam.distThres, mparam.kNN, 1 ) ;
        %catch
        %    pause;
        %end
    end
    %simdot
    
    if bVerbose, fprintf('   %f secs elapsed for matching %d-%d of type %d features\n', toc, length(range1), length(range2), i );   end
        
    if mparam.bMatchDistribution == 2 % equal distribution for each feat type
        nMaxMatchForThis = ceil( mparam.nMaxMatch / nFeatTypeToMatch );    
        [ temp tmpMatchIdx ] = sort(simdot,'descend');
    
        if mparam.bFilterMatch && size(tmpInitialMatch,2) > 0
            % eliminate redundant matches based on the position
            % accumulate features accoring to the rank of sqdist (ascending)
            sel_idx = tmpMatchIdx(1);
            XY_sel_view1 = viewInfo(view1).frame(1:2, range1(tmpInitialMatch(1,tmpMatchIdx(1))));
            XY_sel_view2 = viewInfo(view2).frame(1:2, range2(tmpInitialMatch(2,tmpMatchIdx(1))));
            nSelMatch = 1;
            for  iter_cand = 1:size(tmpInitialMatch,2)
                if nSelMatch >= nMaxMatchForThis, break; end
                XY_view1 = viewInfo(view1).frame(1:2, range1(tmpInitialMatch(1,tmpMatchIdx(iter_cand))));
                XY_view2 = viewInfo(view2).frame(1:2, range2(tmpInitialMatch(2,tmpMatchIdx(iter_cand))));
                if ~bSelf
                    %finds the points within the search radius
                    ridx_view1=BruteSearchMex(XY_sel_view1,XY_view1,'r',mparam.redundancy_thres);
                    ridx_view2=BruteSearchMex(XY_sel_view2,XY_view2,'r',mparam.redundancy_thres);
                else
                    %finds the points within the search radius
                    ridx_view1=BruteSearchMex([XY_sel_view1 XY_sel_view2],XY_view1,'r',mparam.redundancy_thres);
                    ridx_view2=BruteSearchMex([XY_sel_view2 XY_sel_view1],XY_view2,'r',mparam.redundancy_thres);
                end
                equi_ridx = intersect(ridx_view1, ridx_view2);
                if isempty(equi_ridx) % if there's no equi matches
                    nSelMatch = nSelMatch + 1;
                    % insert it into the selected xy list
                    sel_idx(nSelMatch) = tmpMatchIdx(iter_cand);
                    XY_sel_view1(:,nSelMatch) = XY_view1;
                    XY_sel_view2(:,nSelMatch) = XY_view2;
                end            
            end
            if bVerbose, fprintf('   %d candidates selected avoiding equivalent matches (%d pixels)\n', length(sel_idx), mparam.redundancy_thres); end
            sel_idx = sort(sel_idx); % re-aline ascending order
            tmpInitialMatch = tmpInitialMatch(:, sel_idx );
            simdot = simdot( sel_idx );
            
        elseif size(tmpInitialMatch,2) > nMaxMatchForThis
            % select the best nMaxMatch
            if bVerbose, fprintf('   delete %d candidates due to max num of match (%d), max dist:%.2f \n', size(tmpInitialMatch,2)-nMaxMatchForThis, nMaxMatchForThis, temp(nMaxMatchForThis)); end
            del_idx = tmpMatchIdx((nMaxMatchForThis+1):end);
            tmpInitialMatch(:, del_idx ) = [];
            simdot( del_idx ) = [];            
        end            
    end % loop end of equal max distribution
    
    nMatchOfFeatType(i) = size(tmpInitialMatch,2);
    %fprintf('->> %d valid matches for type %d features\n', nMatchOfFeatType(i), i );
    % stack the new matches
    if nMatchOfFeatType(i) > 0
        % initialMatch = [ feat1 feat2 feattype ; ... ]
        initialMatch = [ initialMatch; range1( tmpInitialMatch(1,:) )', range2( tmpInitialMatch(2,:) )', ones(nMatchOfFeatType(i),1)*i ];
        simInitialMatch = [ simInitialMatch simdot ]; 
    end
        
end

nInitialMatches = sum(nMatchOfFeatType);

% When MatchDistribution type == 1, a post-proposs goes on...
% caution! from now on, the match matrix is represented by initialMatch( nInitialMatch, 2 ) 
if mparam.bMatchDistribution == 1 
    [ temp tmpMatchIdx ] = sort(simInitialMatch,'descend');
    if mparam.bFilterMatch && nInitialMatches > 0
        % eliminate euivalent matches based on the position
        % accumulate features accoring to the rank of sqdist (ascending)
        sel_idx(1) = tmpMatchIdx(1);
        XY_sel_view1(1,:) = viewInfo(view1).frame(1:2, initialMatch(tmpMatchIdx(1),1));
        XY_sel_view2(1,:) = viewInfo(view2).frame(1:2, initialMatch(tmpMatchIdx(1),2));
        nSelMatch = 1;
        for  iter_cand = 2:nInitialMatches
            if nSelMatch >= mparam.nMaxMatch, break; end
            XY_view1 = viewInfo(view1).frame(1:2, initialMatch(tmpMatchIdx(iter_cand),1));
            XY_view2 = viewInfo(view2).frame(1:2, initialMatch(tmpMatchIdx(iter_cand),2));
            if ~bSelf
                %finds the points within the search radius
                ridx_view1=BruteSearchMex(XY_sel_view1,XY_view1,'r',mparam.redundancy_thres);
                ridx_view2=BruteSearchMex(XY_sel_view2,XY_view2,'r',mparam.redundancy_thres);
            else
                %finds the points within the search radius
                ridx_view1=BruteSearchMex([XY_sel_view1 XY_sel_view2],XY_view1,'r',mparam.redundancy_thres);
                ridx_view2=BruteSearchMex([XY_sel_view2 XY_sel_view1],XY_view2,'r',mparam.redundancy_thres);
            end
            equi_ridx = intersect(ridx_view1, ridx_view2);
            if isempty(equi_ridx) % if there's no equi matches
                nSelMatch = nSelMatch + 1;
                % insert it into the selected xy list
                sel_idx(nSelMatch) = tmpMatchIdx(iter_cand);
                XY_sel_view1(:,nSelMatch) = XY_view1;
                XY_sel_view2(:,nSelMatch) = XY_view2;
            end            
        end
        if bVerbose, fprintf('   %d candidates selected avoiding equivalent matches (%d pixels)\n', length(sel_idx), mparam.redundancy_thres); end
        sel_idx = sort(sel_idx); % re-aline ascending order
        initialMatch = initialMatch(sel_idx,:);
        simInitialMatch = simInitialMatch( sel_idx );
        nInitialMatches = nSelMatch;
    else 
        if nInitialMatches > mparam.nMaxMatch
            % select the best nMaxMatch 
            if bVerbose, fprintf('   %d matches are eliminated due to max num of match (%d), max dist:%f \n',...
                nInitialMatches-mparam.nMaxMatch, mparam.nMaxMatch, temp(mparam.nMaxMatch)); end
            initialMatch(tmpMatchIdx((mparam.nMaxMatch+1):end),:) = [];
            simInitialMatch(tmpMatchIdx((mparam.nMaxMatch+1):end)) = [];
            nInitialMatches = mparam.nMaxMatch;
        end
    end
    
    
end
%fprintf('>>>> %d total valid matches\n', nInitialMatches );

for k=1:maxFeatType
    nMatch = sum( initialMatch(:,3) == k );
    if nMatch > 0 && bVerbose,  fprintf('     %4d matches from feat type %d\n', nMatch, k ); end
end

%% construct match data
%disp('-- Constructing the initial match data');
matchInfo.match = initialMatch(:,1:2)';
matchInfo.dist = max(simInitialMatch) - simInitialMatch;
matchInfo.sim = simInitialMatch;

if bShow
    showInitialMatches2;
end