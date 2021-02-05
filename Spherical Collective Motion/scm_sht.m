%
% Author: Dhananjay Bhaskar <dhananjay_bhaskar@brown.edu>
% Description: Compute coefficients for SHT to quantify spatiotemporal
% patterns in collective behavior on a sphere
%

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
% Intel IPP Library 7.0 (Realistic Rendering and 3D Data Processing):
%       http://scc.ustc.edu.cn/zlsc/sugon/intel/ipp/ipp_manual/
% Note: The Small Matrices Processing and Realistic Rendering domains
%       containing SHT functions were removed from the Intel IPP version 9.0
%

close all; clear all;

addpath(genpath('Spherical-Harmonic-Transform'));

parent_folder = "SCM_sim_test_5";

% params
N = 10;
maxTime = 400;

coeff = zeros(maxTime, N+1, 2*N+1);

% pick colors for m values
mcolors = distinguishable_colors(2*N+1, 'w');

% create empty figure with legend 
fig1 = figure('Position', [50 50 600 1200]);
hold on;

for l = 0:N
    subplot(N+1,1,l+1)
    xlim([0, 400])
end

% create legend for m values
cnt = 0;
for m = -N:N
    cnt = cnt + 1;
    annotation('textbox', [0.93, 0.92 - cnt*0.015, 0.09, 0.01], 'string', strcat("m = ", num2str(m)), ... 
               'Color', mcolors(cnt,:), 'EdgeColor', 'none')
end
hold off;

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
    % (2l+1) SH functions for each order l (an SH band): -l,...,l
    % (N+1)^2 SH functions for all orders up to N
    [F_N, Y_N] = leastSquaresSHT(N, rho, dirs, 'real');
    
    cnt = 0;
    for l = 0:N
        
        figure(fig1)
        hs = subplot(N+1, 1, l+1);
        cla(hs);
        
        mcount = 0;
        for m = -l:l
            
            cnt = cnt + 1;
            mcount = mcount + 1;
            
            coeff(tp+1,l+1,cnt) = F_N(cnt);
            
            plot(0:tp, coeff(1:(tp+1),l+1,cnt), 'Color', mcolors(mcount,:));
            hold on;
            
        end
        
        hold off;
        xlim([0,400])
        ylabel(strcat("l = ", num2str(l)));
        
    end
    xlabel("Time")
    
    drawnow
    
end

