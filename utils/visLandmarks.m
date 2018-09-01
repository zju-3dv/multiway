function visLandmarks(imgs,W)
N = ceil(sqrt(length(imgs)));
c = mat2gray(sum(W(1:2,:),1));
markersize = 50;
for i = 1:N
    for j = 1:N
        ID = (i-1)*N + j;
        if ID <= length(imgs)
            subplot('position',[(i-1)/N (N-j)/N 1/N 1/N]);
            imshow(imgs{ID}); hold on
            scatter(W(2*ID-1,:),W(2*ID,:),markersize,c,'filled');
            colormap jet;
        end
    end
end