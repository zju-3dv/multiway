function displayMatching(cdata,varargin)

delta = 1;
showFeature = false;
drawEllipse = false;
grayImage = false;
blackBG = true;
linewidth = 1;
color = {'r','g','b','c','y','m'};

ivargin = 1;
while ivargin <= length(varargin)
    switch lower(varargin{ivargin})
        case 'delta'
            ivargin = ivargin + 1;
            delta = varargin{ivargin};
        case 'showfeatures'
            showFeature = true;
        case 'drawellipse'
            drawEllipse = true;
        case 'grayimage'
            grayImage = true;
        case 'linewidth'
            ivargin = ivargin + 1;
            linewidth = varargin{ivargin};
        case 'color'
            ivargin = ivargin + 1;
            color = varargin{ivargin};
        case 'blackbg'
            blackBG = true;
        otherwise
            disp('unknown options, igonored!');
    end
    ivargin = ivargin + 1;
end

for i = 1:2
    load(cdata.viewPath{i});
    view(i) = viewInfo;
end

if isfield(cdata,'X')
    X = cdata.X;
else
    X = [];
end

if isfield(cdata,'GT')
    GT = cdata.GT;
else
    GT = X;
end

%% show features
if showFeature
    for i = 1:2
        subplot(1,2,i);
        imshow(view(i).img); hold on
        visualizeFeatures(view(i).frame,'style','frame');
    end
    return
end

%%
% figure;

imgInput = appendimages(view(1).img,view(2).img);
imgInput = double(imgInput)./255;
if grayImage && ~ismatrix(imgInput)
    imgInput = rgb2gray(imgInput);
end
imshow(imgInput); hold on;
iptsetpref('ImshowBorder','tight');

% draw false matches
curMatchList = cell2mat({cdata.matchInfo(:).match }');
idxFeat1 = curMatchList(1,:);
idxFeat2 = curMatchList(2,:);
feat1 = view(1).frame(:,idxFeat1)';
feat2 = view(2).frame(:,idxFeat2)';
feat2(:,1) = feat2(:,1) + size(view(1).img,2);

% draw true matches
for i = 1:delta:length(X)
    if X(i)
        if GT(i) == 1
            col1 = color{mod(ceil(i/delta),length(color))+1};
            col2 = color{mod(ceil(i/delta),length(color))+1};
        else
            col1 = 'k';
            col2 = 'k';
        end
        
        if blackBG
            plot([ feat1(i,1), feat2(i,1) ]...
                ,[ feat1(i,2), feat2(i,2) ],...
                '-+','LineWidth',linewidth+1,'MarkerSize',10,...
                'color', 'k');
        end
        plot([ feat1(i,1), feat2(i,1) ]...
            ,[ feat1(i,2), feat2(i,2) ],...
            '-+','LineWidth',linewidth,'MarkerSize',10,...
            'color', col1);
        
        if drawEllipse
            drawellipse2( feat1(i,1:5), 1, 'k',5);
            drawellipse2( feat1(i,1:5), 1, col2,4);
            drawellipse2( feat2(i,1:5) ,1, 'k',5);
            drawellipse2( feat2(i,1:5) ,1, col2,4);
        end
    end
end

hold off