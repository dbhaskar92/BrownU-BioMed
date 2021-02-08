%
% Author: Dhananjay Bhaskar <dhananjay_bhaskar@brown.edu>
% Description: Compute coefficients for SHT to quantify spatiotemporal
% patterns in collective behavior on a sphere
%
% References:
%
% MATLAB library by Archontis Politis
%       https://www.mathworks.com/matlabcentral/fileexchange/43856-real-complex-spherical-harmonic-transform-gaunt-coefficients-and-rotations
%       http://research.spa.aalto.fi/projects/sht-lib/sht.html
%
% libsharp:
%       https://www.science.gov/topicpages/s/spherical+harmonic+transform.html
% 
% Interpolation on a sphere:
%       https://web.maths.unsw.edu.au/~rsw/Sphere/#InterpN
%
% Intel IPP Library 7.0 (Realistic Rendering and 3D Data Processing):
%       http://scc.ustc.edu.cn/zlsc/sugon/intel/ipp/ipp_manual/
% Note: The Small Matrices Processing and Realistic Rendering domains
%       containing SHT functions were removed from the Intel IPP version 9.0
%
% FFMPEG commands:
%       ffmpeg -r 10 -f image2 -pattern_type glob -i "*?png" -vcodec libx264 -crf 20 -pix_fmt yuv420p 
%                   -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" output.mp4
%       ffmpeg -r 10 -f image2 -pattern_type glob -i "*?png" -vcodec libx264 -crf 20 -pix_fmt yuv420p 
%                   -filter_complex "[0]pad=w=20+iw:h=20+ih:x=10:y=10:color=black" output.mp4
%       ffmpeg -i large.mp4 -vf "movie=overlay.mp4,scale=350: -1 [inner]; [in][inner] overlay = 800:1210 [out]" output.mp4
%       ffmpeg -i input1.mp4 -i input2.mp4 -filter_complex '[0:v]pad=iw*2:ih[int];[int][1:v]overlay=W/2:0[vid]' 
%                   -map [vid] -c:v libx264 -crf 23 -preset veryfast output.mp4
%

close all; clear all;

addpath(genpath('Spherical-Harmonic-Transform'));

parent_folder = "SCM_sim_test_5";

% params
N = 10;
maxTime = 400;

% output directories
SHT_weights_plot = strcat(parent_folder, filesep, "sht_weights");
SHT_func_approx = strcat(parent_folder, filesep, "sht_func_est");
SHT_weights_folder = convertStringsToChars(SHT_weights_plot);
SHT_func_folder = convertStringsToChars(SHT_func_approx);

if ~exist(SHT_weights_plot, 'dir')
    mkdir(SHT_weights_folder)
end

if ~exist(SHT_func_approx, 'dir')
    mkdir(SHT_func_folder)
end

coeff = zeros(maxTime, N+1, 2*N+1);

% pick colors for m values
mcolors = distinguishable_colors(2*N+1, 'w');
mcolors = circshift(mcolors, N, 1);

% create empty figure with legend 
fig1 = figure('Position', [50 50 600 1200]);
hold on;

for l = 0:N
    subplot(N+1,1,l+1)
    xlim([0, 400])
    ytickformat('%.2f')
end

% create legend for m values
cnt = 0;
for m = -N:N
    cnt = cnt + 1;
    annotation('textbox', [0.93, 0.92 - cnt*0.015, 0.09, 0.01], 'string', strcat("m = ", num2str(m)), ... 
               'Color', mcolors(cnt,:), 'EdgeColor', 'none')
end
hold off;

fig2 = figure('Position',[750 50 300 300]);
set(gca, 'visible', 'off')
xlim([-1.8, 1.8])
ylim([-1.8, 1.8])
zlim([-1.8, 1.8])

for tp = 0:(maxTime-1)
    
    datfile = strcat(parent_folder, filesep, "output_mat", filesep, num2str(tp, '%03.f'), ".mat");

    load(datfile);

    lats = feat_vec(:,1);       % theta, [-pi/2, pi/2]
    lons = feat_vec(:,2);       % phi, [-pi, pi]
    rho = feat_vec(:,3);
    
    % dirs = [azimuth1 inclination1; ...; azimuthK inclinationK]
    % inclination = pi/2 - elevation
    dirs = [lons, pi/2 - lats];
    
    % F_N:  (N+1)^2 vector of SH coefficients
    % Y_N: (N+1)^2 X #dirs, matrix of spherical harmonics evaluated at each
    % (azimuth, inclination) direction
    % Order l : 0,1,...,N
    % (2l+1) SH functions for each order (SH band) l: -l,...,l
    % (N+1)^2 SH functions for all orders up to N
    [F_N, Y_N] = leastSquaresSHT(N, rho, dirs, 'real');
    
    cnt = 0;
    for l = 0:N
        
        figure(fig1)
        hs = subplot(N+1, 1, l+1);
        cla(hs);
        
        for m = -l:l
            
            cnt = cnt + 1;
            
            color_idx = 0;
            for l_val = -N:N
                color_idx = color_idx + 1;
                if m == l_val
                    break
                end
            end
            
            coeff(tp+1,l+1,cnt) = F_N(cnt);
            
            plot(0:tp, coeff(1:(tp+1),l+1,cnt), 'Color', mcolors(color_idx,:));
            hold on;
            
        end
        
        hold off;
        xlim([0,400])
        ylabel(strcat("l = ", num2str(l)));
        ytickformat('%.2f')
        
    end
    xlabel("Time")
    print(strcat(SHT_weights_plot, filesep, sprintf('%03d',tp), ".png"), '-dpng', '-r0')
    
    delta_azi = 1;
    delta_elev = 1;
    
    figure(fig2)
    clf(fig2);
    hax = axes;
    plotSphFunctionCoeffs(F_N, 'real', delta_azi, delta_elev, 'real', hax)
    
    view([75, 21])
    lightangle(75, 0)
    set(gca, 'visible', 'off')
    axis equal
    
    text(1.75, 0, 0, 'X Axis', 'HorizontalAlignment', 'left', 'FontSize', 6, 'Color', 'r');
    text(0, 1.75, 0, 'Y Axis', 'HorizontalAlignment', 'left', 'FontSize', 6, 'Color', 'g');
    text(0, 0, 1.75, 'Z Axis', 'HorizontalAlignment', 'left', 'FontSize', 6, 'Color', 'b');
    
    xlim([-1.75, 1.75])
    ylim([-1.75, 1.75])
    zlim([-1.75, 1.75])
    
    fig2.PaperPositionMode = 'auto';
    print(strcat(SHT_func_approx, filesep, sprintf('%03d',tp), ".png"), '-dpng', '-r0')
    
    drawnow
    
end

save(strcat(parent_folder, filesep, 'SHT_weights.mat'), 'coeff', 'N', 'maxTime', 'parent_folder');
