function acc = evalPMatch(datapath,pMatch,th,mode)
% mode = 1 for PCK, mode = 2 for Precision
if nargin < 3
    mode = 1;
end

viewPath1=[datapath,'/',pMatch.filename{1}];
if isunix()
    viewPath1=strrep(viewPath1,'\','/');
end
view1 = load(viewPath1,'frame','filename');
load(strrep(sprintf('%s/%s.annot.mat',datapath,view1.filename),'\','/'));
key1 = keypoints;
% valid1 = IsInBox(key1,[1 1 size(view1.img,2) size(view1.img,1)]);
% valid1 = valid1(:) & visibility(:);

viewPath2=[datapath,'/',pMatch.filename{2}];
if isunix()
    viewPath2=strrep(viewPath2,'\','/');
end
view2 = load(viewPath2,'frame','filename','img');
load(strrep(sprintf('%s/%s.annot.mat',datapath,view2.filename),'\','/'));
key2 = keypoints;
% valid2 = IsInBox(key2,[1 1 size(view2.img,2) size(view2.img,1)]);
% valid2 = valid2(:) & visibility(:);
% visFlag = valid1 & valid2;
% key1 = key1(:,visFlag);
% key2 = key2(:,visFlag);

idx1 = pMatch.matchInfo.match(1,pMatch.X);
idx2 = pMatch.matchInfo.match(2,pMatch.X);
feat1 = double(view1.frame(1:2,idx1));
feat2 = double(view2.frame(1:2,idx2));

if mode == 2
    tmp = key1;
    key1 = feat1;
    feat1 = tmp;
    tmp = key2;
    key2 = feat2;
    feat2 = tmp;
end

if length(feat1) < 3
    acc = zeros(size(th));
    return
end

F = scatteredInterpolant(feat1(1,:)',feat1(2,:)',feat2(1,:)','natural','nearest');
key2p(1,:) = F(key1(1,:)',key1(2,:)')';
F = scatteredInterpolant(feat1(1,:)',feat1(2,:)',feat2(2,:)','natural','nearest');
key2p(2,:) = F(key1(1,:)',key1(2,:)')';

dist = sqrt(sum((key2-key2p).^2,1));

for i = 1:length(th)
    nTruth(i) = sum(dist<(th(i)*max(size(view2.img))));
end

acc = nTruth / size(key1,2);




