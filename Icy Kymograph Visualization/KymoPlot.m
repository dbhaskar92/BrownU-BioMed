%
% Compute distance travelled and speed of leader cells identified in Icy
% Author: Dhananjay Bhaskar <dhananjay_bhaskar@brown.edu>
% Last Modified: Jun 12, 2018
%

function KymoPlot

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
    
    outfilename = strcat('result', datestr(now,'yyyymmddHHMMSS'), '.csv');
    if size(subfolders, 2) > 0 && ~exist(outfilename, 'file')
        fid = fopen(outfilename, 'wt');
        fprintf(fid, 'Well,CellID,Frame,LineNum,Slope,Dist(microns),Duration(hrs),Speed(microns_per_hr),x1,x2,y1,y2\n');
    else
        return;
    end
    
    for cc = 1 : size(subfolders, 2)
       
        parent_folder = subfolders{cc};
        fprintf('\nWorking Directory: %s\n', parent_folder);
        elements = strsplit(parent_folder, '_');
        
        well_str = char(string(elements(1)));
        well_num = str2double(well_str(2:end));
        
        leader_str = char(string(elements(2)));
        leader_id = str2double(leader_str(2:end));

        x_tif_file = dir(strcat(pwd, filesep, parent_folder, filesep, '*_x.tif'));

        if isempty(x_tif_file)
            disp('Error: TIF file for X not available');
            return;
        end

        x_kymo = imread(strcat(x_tif_file.folder, filesep, x_tif_file.name));
        
        xml_file = dir(strcat(pwd, filesep, parent_folder, filesep, '*t*.xml'));
        
        if isempty(xml_file)
            disp('Error: XML file not available');
            return;
        end
        
        elements = strsplit(xml_file.name, 't');
        elements = strsplit(string(elements(2)), '.');
        time_frame = str2double(char(elements(1)));

        endpt_array = {};

        cnt = 1;
        user_done = false;

        % Ask user to draw lines repeatedly
        while user_done == false

            msg = strcat('Input Line:', {' '}, num2str(cnt));
            fprintf('%s\n', msg{:});
            
            close all;
            figure

            colormap gray

            imagesc(x_kymo);
            colorbar();

            hold on

            h = imline;
            setColor(h, 'green');

            endpts = wait(h);
            endpt_array{end+1} = endpts;

            confirm = input('Are you done? y for Yes, n for No\n', 's');

            delete(h);

            if strcmp(confirm, 'y')
                user_done = true;
            end

            cnt = cnt + 1;

        end
    
        close all;

        % Store figure with lines drawn
        fig = figure;
        colormap gray
        set(fig, 'Visible', 'off');

        title('Kymograph (Horizontal Direction)');
        fig_name = strcat(parent_folder, '_Slope.png');

        imagesc(x_kymo);
        hold on
        colorbar();
    
        for i = 1 : size(endpt_array, 2)

            pts = endpt_array{i};
            xcoords = pts(:, 1);
            ycoords = pts(:, 2);
            plot(xcoords, ycoords, 'b-', 'linewidth', 1.2);

            x_c = mean(pts(:, 1));
            y_c = mean(pts(:, 2));
            anno_str = strcat('Line', {' '}, num2str(i));
            anno = text(x_c, y_c, anno_str);
            anno.Color = [0 1 0];
            anno.FontSize = 8;

            slope = abs((ycoords(2) - ycoords(1))/(xcoords(2)-xcoords(1)));
            dist_px = abs(xcoords(2) - xcoords(1));
            time_frames = abs(ycoords(2) - ycoords(1));

            [dist_um] = px2um(dist_px);
            [time_hrs] = frame2hrs(time_frames);
            speed = dist_um/time_hrs;

            disp('')
            disp(strcat('Line:', {' '}, num2str(i)))
            disp(strcat('|Slope|:', {' '}, num2str(slope)))
            disp(strcat('Speed:', {' '}, num2str(speed), ' um/hrs'))
            disp(strcat('Distance:', {' '}, num2str(dist_um), ' um'))
            disp(strcat('Time:', {' '}, num2str(time_hrs), ' hrs'))
            
            fprintf(fid, '%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f\n', ...
            well_num, leader_id, time_frame, i, slope, dist_um, time_hrs, speed, xcoords(1), xcoords(2), ycoords(1), ycoords(2));

        end

        xlabel('Segment Length ($\mu m$)', 'Interpreter', 'latex');
        ylabel('Time (hrs)', 'Interpreter', 'latex');

        % Display X axis labels in microns
        cur_xt = xticks;
        new_xt_min = min(ceil(px2um(cur_xt)));
        new_xt_max = max(floor(px2um(cur_xt)));
        new_xt = linspace(new_xt_min, new_xt_max, 5);
        xticks(um2px(new_xt));
        xticklabels(new_xt);

        % Display Y axis labels in hours
        cur_yt = yticks;
        new_yt = frame2hrs(cur_yt);
        yticks(hrs2frame(new_yt));
        yticklabels(new_yt);

        set(fig, 'PaperUnits', 'centimeters');
        set(fig, 'PaperPosition', [0 0 7.5 7.5])
        saveas(fig, strcat(parent_folder, filesep, fig_name), 'png');
        
    end
    
    fclose(fid);

end

%
% Convert px to um (microns)
%
function [um_dist_list] = px2um(px_dist_list)

    um_dist_list = zeros(1, size(px_dist_list, 2));
    
    for j = 1 : size(px_dist_list, 2)
        um_dist_list(j) = px_dist_list(j) * 0.65;
    end
    
end

%
% Convert um (microns) to px (pixels)
%
function [px_dist_list] = um2px(um_dist_list)

    px_dist_list = zeros(1, size(um_dist_list, 2));
    
    for j = 1 : size(um_dist_list, 2)
        px_dist_list(j) = um_dist_list(j)/0.65;
    end
    
end

%
% Convert number of frames to hours
%
function [time_elapsed_list] = frame2hrs(num_frames_list)

    time_elapsed_list = zeros(1, size(num_frames_list, 2));
    
    for j = 1 : size(num_frames_list, 2)
        time_elapsed_list(j) = num_frames_list(j) * 0.25;
    end
    
end

%
% Convert hours to number of frames
%
function [num_frames_list] = hrs2frame(time_elapsed_list)

    num_frames_list = zeros(1, size(time_elapsed_list, 2));
    
    for j = 1 : size(time_elapsed_list, 2)
        num_frames_list(j) = time_elapsed_list(j)/0.25;
    end
    
end