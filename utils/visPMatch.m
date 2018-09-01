function visPMatch(datapath,pMatch,mode,varargin)

if isempty(pMatch.X)
    return
end

if nargin < 2
    mode = 1;
end

linewidth = 1;
markersize = 3;
color = {'y'};
bg = false;
th = 0.1;
direction = 'h';

i = 1;
while i <= length(varargin)
    switch lower(varargin{i})
        case 'linewidth'
            i = i+1;
            linewidth = varargin{i};
        case 'markersize'
            i = i+1;
            markersize = varargin{i};
        case 'color'
            i = i+1;
            color = varargin{i};
        case 'th'
            i = i+1;
            th = varargin{i};
        case 'bg'
            i = i+1;
            bg = true;
        case 'direction'
            i = i+1;
            direction = varargin{i};
        otherwise
            fprintf('Unknown option %s is igonored !!!\n',varargin{i});
    end
    i = i+1;
end
            
view1 = load([datapath,'/',pMatch.filename{1}],'img','frame','filename');
view2 = load([datapath,'/',pMatch.filename{2}],'img','frame','filename');
idx1 = pMatch.matchInfo.match(1,pMatch.X);
idx2 = pMatch.matchInfo.match(2,pMatch.X);
feat1 = double(view1.frame(1:2,idx1));
feat2 = double(view2.frame(1:2,idx2));

[imgInput,margin] = appendimages(view1.img,view2.img,direction,10);
% imgInput = mat2gray(imgInput);
imshow(imgInput); hold on;
iptsetpref('ImshowBorder','tight');

feat2v = feat2;
switch direction
    case 'h'
        feat2v(1,:) = feat2v(1,:) + size(view1.img,2) + margin;
    case 'v'
        feat2v(2,:) = feat2v(2,:) + size(view1.img,1) + margin;
end

switch mode
    case 1 % lines
        for i = 1:size(feat1,2)
            if bg
            plot([ feat1(1,i), feat2v(1,i) ]...
                ,[ feat1(2,i), feat2v(2,i) ],...
                '-o','LineWidth',linewidth+1,'MarkerSize',markersize+1,...
                'color', 'k', 'MarkerFaceColor', 'k');
            end
            plot([ feat1(1,i), feat2v(1,i) ]...
                ,[ feat1(2,i), feat2v(2,i) ],...
                '-o','LineWidth',linewidth,'MarkerSize',markersize,...
                'color',color{mod(i,length(color))+1},'MarkerFaceColor',color{mod(i,length(color))+1});
        end
        
    case 2 % color dots
%         c = linspace(0,1,size(feat1,2));
        c = mat2gray(sum(feat1,1)); 
        scatter(feat1(1,:),feat1(2,:),markersize,c,'filled');
        scatter(feat2v(1,:),feat2v(2,:),markersize,c,'filled');
%         for i = 1:length(feat1)
%             text(feat1(1,i),feat1(2,i),num2str(i));
%             text(feat2v(1,i),feat2v(2,i),num2str(i));
%         end
    case 3 % lines with true or false
        if exist(sprintf('%s/%s.mat',datapath,view1.filename(1:end-4)),'file')
            load(sprintf('%s/%s.mat',datapath,view1.filename(1:end-4)));
            key1 = pts_coord;
            load(sprintf('%s/%s.mat',datapath,view2.filename(1:end-4)));
            key2 = pts_coord;
        else
            load(sprintf('%s/%s.pts.mat',datapath,view1.filename));
            key1 = d;
            load(sprintf('%s/%s.pts.mat',datapath,view2.filename));
            key2 = d;
        end
        
        F = scatteredInterpolant(key1(1,:)',key1(2,:)',key2(1,:)','natural','nearest');
        feat2p(1,:) = F(feat1(1,:)',feat1(2,:)')';
        F = scatteredInterpolant(key1(1,:)',key1(2,:)',key2(2,:)','natural','nearest');
        feat2p(2,:) = F(feat1(1,:)',feat1(2,:)')';
        dist = sqrt(sum((feat2-feat2p).^2,1));
        idxTrue = dist < th * max(size(view2.img));
        
        for i = find(idxTrue)
            if bg
                plot([ feat1(1,i), feat2v(1,i) ]...
                                ,[ feat1(2,i), feat2v(2,i) ],...
                                '-+','LineWidth',1,'MarkerSize',10,...
                                'color', 'k');
            end
            plot([ feat1(1,i), feat2v(1,i) ]...
                ,[ feat1(2,i), feat2v(2,i) ],...
                '-+','LineWidth',1,'MarkerSize',10,...
                'color', 'b');
        end
        
        for i = find(~idxTrue)
            plot([ feat1(1,i), feat2v(1,i) ]...
                ,[ feat1(2,i), feat2v(2,i) ],...
                '-+','LineWidth',2,'MarkerSize',10,...
                'color', 'r');
        end
end












