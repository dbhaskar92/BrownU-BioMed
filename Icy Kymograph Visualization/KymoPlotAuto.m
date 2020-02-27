%
% Plot kymographs quantifying collective migration and cluster aggregation
% Author: Dhananjay Bhaskar <dhananjay_bhaskar@brown.edu>
% Last Modified: Jun 7, 2018
%

function KymoPlotAuto

    % Get directories in current folder
    cwd_contents = dir(pwd);
    dir_indx = [cwd_contents.isdir];
    cwd_dirs = cwd_contents(dir_indx);
    subfolders = {};
    
    for cnt = 1 : size(cwd_dirs)
       
        if strcmp(cwd_dirs(cnt).name, '.')
            continue;
        elseif strcmp(cwd_dirs(cnt).name, '..')
            continue;
        else
            subfolders{end+1} = cwd_dirs(cnt).name;
        end
        
    end
    
    % Process each directory
    for i = 1 : size(subfolders, 2)
       
        folder = subfolders{i};
        
        msg = strcat('Processing:', {' '}, num2str(i), '/', num2str(size(subfolders, 2)), ' (', folder, ')');
        fprintf('%s\n', msg{:});
        
        x_tif_file = dir(strcat(pwd, filesep, folder, filesep, '*_x.tif'));
        y1_tif_file = dir(strcat(pwd, filesep, folder, filesep, '*_y1.tif'));
        y2_tif_file = dir(strcat(pwd, filesep, folder, filesep, '*_y2.tif'));

        y_files_exist = 1;
        
        if isempty(y1_tif_file)
            y_files_exist = 0;
        end
        if isempty(y2_tif_file)
            y_files_exist = 0;
        end

        x_kymo = imread(strcat(x_tif_file.folder, filesep, x_tif_file.name));
        num_frames = size(x_kymo, 1);
        
        if y_files_exist == 1
            
            y1_kymo = imread(strcat(y1_tif_file.folder, filesep, y1_tif_file.name));
            y2_kymo = imread(strcat(y2_tif_file.folder, filesep, y2_tif_file.name));
        
            if (num_frames ~= size(y1_kymo, 1) || num_frames ~= size(y2_kymo, 1))
                disp('Error: Number of frames does not match.');
                exit();
            end
            
        end
        
        time_elapsed = frame2hrs(num_frames);
        msg = strcat('Number of frames:', {' '}, num2str(num_frames), ' (', num2str(time_elapsed), ' hrs)');
        fprintf('%s\n', msg{:});
        
        x_segment_length = px2um(size(x_kymo, 2));
        msg = strcat('Length of x segment:', {' '}, num2str(x_segment_length), ' microns');
        fprintf('%s\n', msg{:});
        
        if y_files_exist == 1
            y1_segment_length = px2um(size(y1_kymo, 2));
            msg = strcat('Length of y1 segment:', {' '}, num2str(y1_segment_length), ' microns');
            fprintf('%s\n', msg{:});

            y2_segment_length = px2um(size(y2_kymo, 2));
            msg = strcat('Length of y2 segment:', {' '}, num2str(y2_segment_length), ' microns');
            fprintf('%s\n', msg{:});
        end
        
        img_file = dir(strcat(pwd, filesep, folder, filesep, '*t*.tif'));
        raw_img = imread(strcat(img_file.folder, filesep, img_file.name));
        
        % Parse XML file
        xml_file = dir(strcat(pwd, filesep, folder, filesep, '*t*.xml'));
        xml = parseXML(strcat(xml_file.folder, filesep, xml_file.name));
        
        % Extract (x_0, y_0) and (x_1, y_1) for each line segment
        [l1_X, l1_Y, l2_X, l2_Y, l3_X, l3_Y] = find_lines(xml, y_files_exist);
        
        fig1 = figure;
        set(fig1, 'Visible', 'off');
        title('Kymograph Plots');
        fig1_name = strcat('KymoP_', folder);
        colormap gray
        
        if y_files_exist == 1
            
            subplot(3,1,1)
            imagesc(x_kymo)
            colorbar()
            xlabel('X Index')
            ylabel('Frame')
            
            ax = gca;
            ax.XColor = 'black';
            ax.YColor = 'black';
            ax.YGrid = 'on';
            set(gca, 'fontsize', 5);

            subplot(3,1,2)
            imagesc(y1_kymo)
            colorbar()
            xlabel('Y1 Index')
            ylabel('Frame')
            
            ax = gca;
            ax.XColor = 'black';
            ax.YColor = 'black';
            ax.YGrid = 'on';
            set(gca, 'fontsize', 5);

            subplot(3,1,3)
            imagesc(y2_kymo)
            colorbar()
            xlabel('Y2 Index')
            ylabel('Frame')
            
            ax = gca;
            ax.XColor = 'black';
            ax.YColor = 'black';
            ax.YGrid = 'on';
            set(gca, 'fontsize', 5);
            
        else
            
            imagesc(x_kymo)
            colorbar()
            xlabel('X Index')
            ylabel('Frame')
            
            ax = gca;
            ax.XColor = 'black';
            ax.YColor = 'black';
            ax.YGrid = 'on';
            set(gca, 'fontsize', 5);
            
        end
        
        set(fig1, 'PaperUnits', 'centimeters');
        set(fig1, 'PaperPosition', [0 0 7.5 7.5])
        
        saveas(fig1, strcat(pwd, filesep, folder, filesep, fig1_name), 'epsc');
        print('-dtiff', '-r1000', strcat(pwd, filesep, folder, filesep, fig1_name)) 

        fig2 = figure;
        set(fig2, 'Visible', 'off');
        title('Kymograph Lines');
        fig2_name = strcat('KymoM_', folder);
        
        imshow(imadjust(raw_img))
        hold on

        plot(l1_X, l1_Y, 'g-', 'linewidth', 1.1)
        plot(l2_X, l2_Y, 'g-', 'linewidth', 1.1)
        plot(l3_X, l3_Y, 'g-', 'linewidth', 1.1)

        ax = gca;
        ax.XColor = 'black';
        ax.YColor = 'black';
        ax.YGrid = 'on';
        set(gca, 'fontsize', 1);
        set(gca, 'XMinorTick', 'off', 'YMinorTick', 'off')

        set(fig2, 'PaperUnits', 'centimeters');
        set(fig2, 'PaperPosition', [0 0 7.5 7.5])
        saveas(fig2, strcat(pwd, filesep, folder, filesep, fig2_name), 'epsc');
        print('-dtiff', '-r1000', strcat(pwd, filesep, folder, filesep, fig2_name))
        
    end
    
end

%
% Parse XML file to extract (x,y) coordinates of line segment endpoints
% Reference:
% xml.Children(6).Name = rois
% xml.Children(6).Children(idx).Name = 'roi', idx = 2, 4, 6
% xml.Children(6).Children(2).Children(28).Name = 'pt1'
% xml.Children(6).Children(2).Children(30).Name = 'pt2'
% xml.Children(6).Children(2).Children(28).Children(2).Name = 'pos_x'
% xml.Children(6).Children(2).Children(28).Children(4).Name = 'pos_y'
%
function [l1_X, l1_Y, l2_X, l2_Y, l3_X, l3_Y] = find_lines(xml, y_files_flag)

    l1_X = [];
    l1_Y = [];
    l2_X = [];
    l2_Y = [];
    l3_X = [];
    l3_Y = [];

    % Find element storing ROI information in XML tree structure
    root_num_children = size(xml.Children, 2);
    roi_data_indx = -1;

    for idx = 1 : root_num_children
        if strcmp(xml.Children(idx).Name, 'rois')
            roi_data_indx = idx;
            break;
        end
    end

    if roi_data_indx == -1
        disp('Error: Could not find ROI data');
        exit();
    end

    roidata_num_children = size(xml.Children(roi_data_indx).Children, 2);
    roi_indices = [];
    
    for idx = 1 : roidata_num_children
        
        if strcmp(xml.Children(roi_data_indx).Children(idx).Name, 'roi')
            
            if size(roi_indices, 2) >= 3
                 disp('Error: Found more than 3 ROIs');
                 exit();
            else    
                roi_indices = [roi_indices idx];
            end
            
        end
        
    end
    
    if (size(roi_indices, 2) < 3 && y_files_flag == 1)
        disp('Error: Did not find 3 ROIs');
        exit();
    end
    
    % For each ROI, find (x,y) coordinates for line segment
    
    num_rois = 3;
    if y_files_flag == 0
        num_rois = 1;
    end
    
    for cnt = 1 : num_rois
        
        roi_idx = roi_indices(cnt);
        roi_num_children = size(xml.Children(roi_data_indx).Children(roi_idx).Children, 2);
        
        pt1_idx = -1;
        pt2_idx = -1;
        
        for idx = 1 : roi_num_children
        
            if strcmp(xml.Children(roi_data_indx).Children(roi_idx).Children(idx).Name, 'pt1')
                pt1_idx = idx;
            elseif strcmp(xml.Children(roi_data_indx).Children(roi_idx).Children(idx).Name, 'pt2') 
                pt2_idx = idx;
            end

        end
        
        if (pt1_idx == -1 || pt2_idx == -1)
            disp('Error: Unable to find line segment endpoints');
            exit();
        end
        
        x1 = xml.Children(roi_data_indx).Children(roi_idx).Children(pt1_idx).Children(2).Children.Data;
        y1 = xml.Children(roi_data_indx).Children(roi_idx).Children(pt1_idx).Children(4).Children.Data;
        x2 = xml.Children(roi_data_indx).Children(roi_idx).Children(pt2_idx).Children(2).Children.Data;
        y2 = xml.Children(roi_data_indx).Children(roi_idx).Children(pt2_idx).Children(4).Children.Data;
        
        if cnt == 1     
            l1_X = [sscanf(x1, '%g,', [1, inf]).' sscanf(x2, '%g,', [1, inf]).'];
            l1_Y = [sscanf(y1, '%g,', [1, inf]).' sscanf(y2, '%g,', [1, inf]).'];
        elseif cnt == 2
            l2_X = [sscanf(x1, '%g,', [1, inf]).' sscanf(x2, '%g,', [1, inf]).'];
            l2_Y = [sscanf(y1, '%g,', [1, inf]).' sscanf(y2, '%g,', [1, inf]).'];
        elseif cnt == 3
            l3_X = [sscanf(x1, '%g,', [1, inf]).' sscanf(x2, '%g,', [1, inf]).'];
            l3_Y = [sscanf(y1, '%g,', [1, inf]).' sscanf(y2, '%g,', [1, inf]).'];
        end
        
    end

end

%
% Convert px to um (microns)
%
function [um_dist] = px2um(px_dist)
    um_dist = px_dist * 0.65;
end

%
% Convert number of frames to hours
%
function [time_elapsed] = frame2hrs(num_frames)
    time_elapsed = num_frames * 0.25;
end