function [modeEnd GT_M GT_T] = make_groundTruth_feat( cdata, GT_M, GT_T )
%% Threshold Parameter
distClose = 10;

%% Monitor Parameter
ScreenResolution = get(0,'ScreenSize');
UserMonitor.Resolution = ScreenResolution(3:4); % Monitor Resolution
UserMonitor.bDualDisplay = 0; % Boolean for Dual Dispaly
%UserMonitor.primaryDisplay = 'left'; % Which one is primary dispaly
UserMonitor.primaryDisplay = 'right'; % Which one is primary dispaly

if strcmp(UserMonitor.primaryDisplay, 'right')
    rectM = [0 0 UserMonitor.Resolution./2];
    rectT = [UserMonitor.Resolution(1)/2 0 UserMonitor.Resolution/2];
    rectMenu = [UserMonitor.Resolution(1)/3 0 UserMonitor.Resolution/4];
else
    rectM = [];%[UserMonitor.Resolution(1) 0 UserMonitor.Resolution(1) UserMonitor.Resolution(2)];
    rectT = [];%[0 0 UserMonitor.Resolution];
    rectMenu = [];
end

%% Generate Figures
h_model = figure('NumberTitle', 'Off', 'Name', 'GT Reference', 'Position', rectM);
h_test = figure('NumberTitle', 'Off', 'Name', 'GT Input', 'Position', rectT);
h_all = figure('NumberTitle', 'Off', 'Name', 'GT Results');
vis_img = appendimages(rgb2gray(cdata.view(1).img),rgb2gray(cdata.view(2).img));
gray_img1 = rgb2gray(cdata.view(1).img);
gray_img2 = rgb2gray(cdata.view(2).img);


% Construct Correspondence data for features
frameM = cdata.view(1).frame;
frameT = cdata.view(2).frame;

bStart = 1;
bContinue = 1;
colorCode = makeColorCode(100);

for i=1:size(GT_M,2)
    idx = nearestneighbour( GT_M(1:2,i), frameM(1:2,:), 'Radius', distClose);
    if ~isempty(idx)
        GT_idxFrame_M(i) = idx(1);
    else
        GT_idxFrame_M(i) = 0;
    end
end
for i=1:size(GT_T,2)
    idx = nearestneighbour( GT_T(1:2,i), frameT(1:2,:), 'Radius', distClose);
    if ~isempty(idx)
        GT_idxFrame_T(i) = idx(1);
    else
        GT_idxFrame_T(i) = 0;
    end
end
                
%% Main Loop
while bContinue % Select Feature in Model image
    %% Show current GT
    figure(h_all);
    imshow(vis_img); hold on;
    for i = 1:size(GT_M,2)
        %frame1 = cdata.view(1).frame(:,TEMP(i,1));
        %frame2 = cdata.view(2).frame(:,TEMP(i,2)); 
        %frame2(1) = frame2(1) + size(cdata.view(1).img,2);
        if i <= size(GT_T,2)
         plot([ GT_M(1,i), GT_T(1,i)+size(cdata.view(1).img,2) ]...
            ,[ GT_M(2,i), GT_T(2,i) ],...
            '-','LineWidth',3,'MarkerSize',10,...
            'color', colorCode(:,i));
        end
        plot(GT_M(1,i), GT_M(2,i), 'o', 'color', colorCode(:,i), 'linestyle', 'none', 'lineWidth', 3, 'MarkerSize', 10);
        if i <= size(GT_T,2)
            plot(GT_T(1,i)+size(cdata.view(1).img,2), GT_T(2,i), 'o', 'color', colorCode(:,i), 'linestyle', 'none', 'lineWidth', 3, 'MarkerSize', 10);
        end
        %visualizeFeatures( feat1, 'style', 'frame', 'colorcode', 'y');
        %visualizeFeatures( feat2, 'style', 'frame', 'colorcode', 'y');
        %clear TEMP
    end
    if bStart
        bStart = 0;
        valMenu = 1;
    else
        valMenu = menu('                            Choose one action                            ', ...
        '                            Input Edition                            ', ...
        '                            Ref. Edition                            ', ...    
        '                            Refresh                            ', ...
        '                            Terminate                           ',...
        '                            Save & Previous                        ',...
        '                            Save & Next                      ');
        
    end
    
    switch valMenu
        %% Case 1: Addition
        case 1
            disp('Left click for GT. Right click for MENU');
            while 1
                figure(h_model); imshow(gray_img1); set(h_model, 'Position', rectM);  hold on;
                % show already selected points
                for i=1:size(GT_M,2)
                    if GT_idxFrame_M(i) > 0
                        visualizeFeatures( frameM(:,GT_idxFrame_M(i)), 'style', 'frame', 'colorcode', colorCode(:,i));
                    end
                    plot(GT_M(1,i), GT_M(2,i), 'o', 'color', colorCode(:,i), 'linestyle', 'none', 'lineWidth', 3, 'MarkerSize', 10);                    
                end
                for i=1:size(GT_M,2)
                    text(GT_M(1,i), GT_M(2,i), num2str(i),...
                        'color','k','HorizontalAlignment', 'center','FontSize',15);
                end
                figure(h_test); imshow(gray_img2); set(h_test, 'Position', rectT); hold on;
                plot(frameT(1,:), frameT(2,:), 'ko', 'linestyle', 'none');
                plot(frameT(1,:), frameT(2,:), 'y*', 'linestyle', 'none');
                for i=1:size(GT_T,2)
                    if GT_idxFrame_T(i) > 0
                        visualizeFeatures( frameT(:,GT_idxFrame_T(i)), 'style', 'frame', 'colorcode', colorCode(:,i));
                    end
                    plot(GT_T(1,i), GT_T(2,i), 'o', 'color', colorCode(:,i), 'linestyle', 'none', 'lineWidth', 3, 'MarkerSize', 10);
                end
                for i=1:size(GT_T,2)
                    text(GT_T(1,i), GT_T(2,i), num2str(i),...
                        'color','k','HorizontalAlignment', 'center','FontSize',15);
                end
                if size(GT_M,2) < 1 
                    disp('No GT in the reference points');
                    break;
                end
                
                %GT_M
                %GT_T
                if size(GT_T,2) > 0
                    idGT = find(isnan(GT_T(1,:)),1,'first');
                    if isempty(idGT)
                        idGT = size(GT_T,2) + 1;
                    end
                    if idGT > size(GT_M,2)
                        idGT = 0;
                    end
                else
                    idGT = 1;
                end
                
                %idGT
                if idGT > 0
                    figure(h_model);                
                    plot(GT_M(1,idGT), GT_M(2,idGT), 'o', 'color', colorCode(:,idGT), 'linestyle', 'none', 'lineWidth', 5, 'MarkerSize', 40);
                    figure(h_test);
                    if size(GT_T,2) >= idGT                        
                        plot(GT_T(1,idGT), GT_T(2,idGT), 'o', 'color', colorCode(:,idGT), 'linestyle', 'none', 'lineWidth', 5, 'MarkerSize', 40);
                    end
                end
                
                figure(h_test);
                [x y button] = ginput(1); % Mouse Input
                if button ~= 1, break; end
                mode = 1; % 1:add 2:erase
                % Find closest gt in Model
                if size(GT_T,2) > 0
                    idx = nearestneighbour( [x; y], GT_T(1:2,:), 'Radius', distClose);
                    if ~isempty(idx)
                        mode = 2;
                        x = GT_T(1,idx(1));
                        y = GT_T(2,idx(1));
                    end
                end
                if mode == 1 && idGT > 0% add
                    % update the GT
                    GT_T(:,idGT) = [x; y];
                    idx = nearestneighbour( [x; y], frameT(1:2,:), 'Radius', distClose);
                    if ~isempty(idx) % update the closest frame
                        GT_idxFrame_T(idGT) = idx(1);                       
                    else
                        GT_idxFrame_T(idGT) = 0;
                    end                   
                elseif mode == 2 % erase
                   GT_T(:,idx(1)) = [nan; nan];
                   GT_idxFrame_T(idx(1)) = 0;
                end
            end
                
        case 2
            disp('Left click for GT. Right click for MENU');
            while 1
                figure(h_model); imshow(gray_img1); set(h_model, 'Position', rectM);  hold on;
                % show already selected points
                for i=1:size(GT_M,2)
                    if GT_idxFrame_M(i) > 0
                        visualizeFeatures( frameM(:,GT_idxFrame_M(i)), 'style', 'frame', 'colorcode', colorCode(:,i));
                    end
                    plot(GT_M(1,i), GT_M(2,i), 'o', 'color', colorCode(:,i), 'linestyle', 'none', 'lineWidth', 3, 'MarkerSize', 10);                    
                end
                for i=1:size(GT_M,2)
                    text(GT_M(1,i), GT_M(2,i), num2str(i),...
                        'color','k','HorizontalAlignment', 'center','FontSize',15);
                end
                plot(frameM(1,:), frameM(2,:), 'ko', 'linestyle', 'none');
                plot(frameM(1,:), frameM(2,:), 'y*', 'linestyle', 'none');
                
                figure(h_model);
                [x y button] = ginput(1); % Mouse Input
                if button ~= 1, break; end
                mode = 1; % 1:add 2:erase
                % Find closest gt in Model
                if size(GT_M,2) > 0
                    idx = nearestneighbour( [x; y], GT_M(1:2,:), 'Radius', distClose);
                    if ~isempty(idx) 
                        mode = 2;
                        x = GT_M(1,idx(1));
                        y = GT_M(2,idx(1));
                    end
                end
                if mode == 1 % add
                   idGT = size(GT_M,2) + 1;
                   idx = nearestneighbour( [x; y], frameM(1:2,:), 'Radius', distClose);
                   if ~isempty(idx)
                        GT_idxFrame_M(idGT) = idx(1);
                        %frame1 = frameM(:,idx);
                        %plot(frame1(1), frame1(2), 'r*');
                        %visualizeFeatures( frame1, 'style', 'frame', 'colorcode', 'r');
                    else
                        GT_idxFrame_M(idGT) = 0;
                    end
                   %plot(x, y, 'ro', 'linestyle', 'none', 'lineWidth', 3, 'MarkerSize', 10);
                   GT_M(:,idGT) = [x; y];
                elseif mode == 2 % erase
                   idGT = idx;
                   %plot(x, y, 'ko', 'linestyle', 'none', 'lineWidth', 3, 'MarkerSize', 10);
                   GT_M(:,idGT) = [];
                   GT_idxFrame_M(idGT) = [];
                end
            end
            
        %% Refresh
        case 3
            figure(h_all); clf;
            figure(h_model); clf;
            figure(h_test); clf;
            
        %% Terminate Program
        case 4 % save and go next
            modeEnd = 3;
            break; % save and go previous
        case 5
            modeEnd = 2;
            break;
        case 6 % terminate
            modeEnd = 1;
            break;    
    end
end

close all; drawnow;

