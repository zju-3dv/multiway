function [feat desc]=loadfeatures_v3(file, detector, thresholdlevel)
% Load local features

outname = [ 'tmp_' detector '_features.out' ]; %% temporary output file buffer name
fpath = fileparts(which(mfilename)); %% current directory path
inname = file;

% Default executable file settings for feature detection 
% by Mikolajczyk et al.'s http://www.featurespace.org/
if strncmp(computer,'PC',2) % MS Windows
  exec_str = ['"' fpath '/detectors/compute_descriptors.exe"'];
elseif strcmp(computer,'GLNX86') % Linux
  exec_str = [fpath '/detectors/compute_descriptors_32bit.ln'];
elseif strcmp(computer,'GLNXA64')% Linux64
  exec_str = [fpath '/detectors/compute_descriptors_64bit.ln'];
else error('This function can run only with MS Windows or Linux');
end

%fprintf('-- detecting %s features from %s\n',detector, file);

desc = [];

switch detector    
        
    case 'mser_b'
        if 1
            %if (nargin<3)
            %   fprintf('No settings specified, so using defaults...\n');
               param.MSER_Ellipse_Scale = 1.0;
               param.MSER_Maximum_Relative_Area = 0.010;
               param.MSER_Minimum_Size_Of_Output_Region = 30;
               param.MSER_Minimum_Margin = 10;
               param.MSER_Use_Relative_Margins = 0;
               param.MSER_Vervose_Output = 0;
            %end
            if thresholdlevel == 2
                param.MSER_Minimum_Size_Of_Output_Region = 15;  % MSER 
                param.MSER_Minimum_Margin = 3;          % MSER  
            else
                param.MSER_Minimum_Size_Of_Output_Region = 15;  % MSER 
                param.MSER_Minimum_Margin = 0; % or 10          % MSER
            end
             
            
            opt = sprintf('-t 2 -es %f -per %f -ms %d -mm %d -i "%s" -o "%s"',...
                param.MSER_Ellipse_Scale, param.MSER_Maximum_Relative_Area,...
                param.MSER_Minimum_Size_Of_Output_Region, param.MSER_Minimum_Margin,...
                inname, outname );
            if strncmp(computer,'PC',2) % MS Windows
              exec_str = ['"' fpath '/detectors/mser.exe"'];
            elseif strcmp(computer,'GLNX86') || strcmp(computer,'GLNXA64')% Linux
              exec_str = [fpath '/detectors/mser.ln'];
            else error('This function can run only with MS Windows or Linux');
            end
        else
            opt = sprintf('-mser -i "%s" -o1 "%s"',inname ,outname );
        end
    
    case 'hesaff_b'
        
        if 1
            if strncmp(computer,'PC',2) % MS Windows
              exec_str = ['"' fpath '/detectors/extract_features_32bit.exe"'];
            elseif strcmp(computer,'GLNX86') % Linux
              exec_str = [fpath '/detectors/extract_features_32bit.ln'];
            elseif strcmp(computer,'GLNXA64')% Linux64
              exec_str = [fpath '/detectors/extract_features_64bit.ln'];
              %exec_str = [fpath '/detectors/detect_points.ln'];  
            else error('This function can run only with MS Windows or Linux');
            end
        end
        
        if thresholdlevel == 2
            param.HESAFF_hesThres = 500;
        else
            %param.HESAFF_hesThres = 10;
            param.HESAFF_hesThres = 100;
        end
        
        opt = sprintf('-hesaff -hesThres %d -i "%s" -o1 "%s"',...
            param.HESAFF_hesThres, ...
            inname, outname );

        

    case 'haraff_b'
        
        if 1
            if strncmp(computer,'PC',2) % MS Windows
              exec_str = ['"' fpath '/detectors/extract_features_32bit.exe"'];
            elseif strcmp(computer,'GLNX86') %Linux
              exec_str = [fpath '/detectors/extract_features_32bit.ln'];
            elseif strcmp(computer,'GLNXA64')% Linux64
              exec_str = [fpath '/detectors/extract_features_64bit.ln'];
              %exec_str = [fpath '/detectors/detect_points.ln'];  
            else error('This function can run only with MS Windows or Linux');
            end
        end

        if thresholdlevel == 2
            param.HARAFF_harThres = 1500;
        else
            %param.HARAFF_harThres = 10;
            param.HARAFF_harThres = 500;
        end

        opt = sprintf('-haraff -harThres %d -i "%s" -o1 "%s"',...
            param.HARAFF_harThres, ...
            inname, outname );

    case 'heslap_b'
        
        if 1
            if strncmp(computer,'PC',2) % MS Windows
              exec_str = ['"' fpath '/detectors/"'];
            elseif strcmp(computer,'GLNX86') %Linux
              exec_str = [fpath '/detectors/h_affine.ln'];
            elseif strcmp(computer,'GLNXA64')% Linux64
              exec_str = [fpath '/detectors/h_affine.ln'];
              %exec_str = [fpath '/detectors/detect_points.ln'];  
            else error('This function can run only with MS Windows or Linux');
            end
        end
        
        if thresholdlevel == 2
            param.HESAFF_hesThres = 300;
        else
            %param.HESAFF_hesThres = 10;
            param.HESAFF_hesThres = 100;
        end
        
        opt = sprintf('-heslap -thres %d -i "%s" -o1 "%s"',...
            param.HESAFF_hesThres, inname, outname );
        
    case 'harlap_b'
        
        if 1
            if strncmp(computer,'PC',2) % MS Windows
              exec_str = ['"' fpath '/detectors/"'];
            elseif strcmp(computer,'GLNX86') %Linux
              exec_str = [fpath '/detectors/h_affine.ln'];
            elseif strcmp(computer,'GLNXA64')% Linux64
              exec_str = [fpath '/detectors/h_affine.ln'];
              %exec_str = [fpath '/detectors/detect_points.ln'];  
            else error('This function can run only with MS Windows or Linux');
            end
        end
        
        if thresholdlevel == 2
            param.HARAFF_harThres = 1500;
        else
            %param.HARAFF_harThres = 10;
            param.HARAFF_harThres = 500;
        end
        
        opt = sprintf('-harlap -thres %d -i "%s" -o1 "%s"',...
            param.HARAFF_harThres, inname, outname );
        
    case 'dog_b' % DOG features
        if 1
            
            if thresholdlevel == 2
                param.DOG_dogThres = 1500;
            else
                %param.HARAFF_harThres = 10;
                param.DOG_dogThres = 300;
            end
            
            opt = sprintf('-dog -thres %d -i "%s" -o1 "%s"',...
                param.DOG_dogThres, inname, outname );
        else
            % Load image
            image = imread(inname);
            if size(image,3) > 1
               image = rgb2gray(image);
            end
            [rows, cols] = size(image); 
            % Convert into PGM imagefile, readable by "keypoints" executable
            f = fopen('tmp.pgm', 'w');
            if f == -1
                error('Could not create file tmp.pgm.');
            end
            fprintf(f, 'P5\n%d\n%d\n255\n', cols, rows);
            fwrite(f, image', 'uint8');
            fclose(f);
            fprintf('-- Detecting SIFT LOG features from %s\n',file);

            opt = sprintf('<tmp.pgm >"%s"', outname);
            if strncmp(computer,'PC',2) % MS Windows
              exec_str = ['"' fpath '/detectors/siftWin32"'];
            elseif strcmp(computer,'GLNX86') % Linux
              exec_str = [fpath '/detectors/sift'];
            elseif strcmp(computer,'GLNXA64')% Linux64
              exec_str = [fpath '/detectors/sift'];  
            else error('This function can run only with MS Windows or Linux');
            end
        end
    
        %if f ~= -1, delete('tmp.pgm');  end
end
        
% Call the binary executable
%[exec_str  ' ' opt ]
result = unix([exec_str  ' ' opt ]);

if result ~= 0
  error('Calling the [ %s ] feature detector failed processing %s.',detector, inname);
end

% Load the output file
fid = fopen(outname, 'r');
if fid==-1
  error('Cannot load results from [ %s ] feature detector processing %s.', detector,inname);
end

if strcmp(detector,'sift')
    [header, count] = fscanf(fid, '%d %d', [1 2]);
    if count ~= 2
        error('Invalid keypoint file beginning.');
    end
    num = header(1);
    len = header(2);
    if len ~= 128
        error('Keypoint descriptor length invalid (should be 128).');
    end

    % Creates the two output matrices (use known size for efficiency)
    %locs = double(zeros(num, 4));
    %descriptors = double(zeros(num, 128));
    feat = zeros(num, 5);
    % Parse tmp.key
    for i = 1:num
        [vector, count] = fscanf(fid, '%f %f %f %f', [1 4]); %row col scale ori
        if count ~= 4
            error('Invalid keypoint file format');
        end
        % convert the params into elliptical representation
        feat(i, 1) = vector(1, 2);
        feat(i, 2) = vector(1, 1);
        feat(i, 3) = 0.05/vector(1, 3)^2; % adequate param: 0.05 & x2 or 0.03 & x1.5
        feat(i, 4) = 0;
        feat(i, 5) = 0.05/vector(1, 3)^2;

        [descrip, count] = fscanf(fid, '%d', [1 len]);
        if (count ~= 128)
            error('Invalid keypoint file value.');
        end
        %Normalize each input vector to unit length
        %descrip = descrip / sqrt(sum(descrip.^2));
        %desc(i, :) = descrip(:);
    end
    
else
    try
        header=fscanf(fid, '%f',1);
        num=fscanf(fid, '%d',1);
        feat = fscanf(fid, '%f', [5, inf]);
        feat=feat';        
        %s = reshape( fdata(c+[1:5*fdata(2)]), [5 fdata(2)] );
    catch
        error('Wrong length of the output file processing %s.',inname);
    end
    
end

fclose(fid);



% %% descriptor extraction
% descriptor = 'sift';
% 
% outname_desc = [ 'tmp_' detector '_' descriptor '.out' ]; %% temporary output file buffer name
% 
% % Default executable file settings for feature detection 
% % by Mikolajczyk et al.'s http://www.featurespace.org/
% if strncmp(computer,'PC',2) % MS Windows
%   exec_str = ['"' fpath '/detectors/compute_descriptors.exe"'];
% elseif strcmp(computer,'GLNX86') % Linux
%   exec_str = [fpath '/detectors/compute_descriptors_32bit.ln'];
% elseif strcmp(computer,'GLNXA64')% Linux64
%   exec_str = [fpath '/detectors/compute_descriptors_64bit.ln'];
% else error('This function can run only with MS Windows or Linux');
% end
% 
% fprintf('-- making %s descriptors from %s\n',descriptor, file);
% 
% opt = sprintf('-sift -noangle -i "%s" -p1 "%s" -o1 "%s"',inname, outname, outname_desc );
% 
% % Call the binary executable
% result = unix([exec_str  ' ' opt ]);
% 
% if result ~= 0
%   error('Calling the [ %s ] feature descriptor failed processing %s.',descriptor, inname);
% end
% 
% % Load the output file
% fid = fopen(outname_desc, 'r');
% if fid==-1
%   error('Cannot load results from [ %s ] feature detector processing %s.', detector,inname);
% end
% 
% try
%     header=fscanf(fid, '%f',1);
%     num=fscanf(fid, '%d',1);
%     desc = fscanf(fid, '%f', [133, inf]);
%     desc = desc(6:end,:)';        
%     %s = reshape( fdata(c+[1:5*fdata(2)]), [5 fdata(2)] );
% catch
%     error('Wrong length of the output file processing %s.',inname);
% end
%     
% fclose(fid);
% 
% desc = desc./repmat(sqrt(sum(desc.^2,2)),[1 128]);

