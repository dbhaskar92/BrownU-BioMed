%
% Description: Quantify image intensity in green channel along 
% user-defined curve
% Author: Dhananjay Bhaskar <dhananjay_bhaskar@brown.edu>
%

close all; clear all; clc;

addpath('bfmatlab');
[imgfile, imgpath] = uigetfile('*.zvi', 'Select a file');
imgdata = bfopen(strcat(imgpath, imgfile));

num_series = size(imgdata, 1);
assert(num_series == 1, 'ERROR: Expected single series acquisition.')

num_planes = size(imgdata{1, 1}, 1);
assert(num_planes == 2, 'ERROR: Expected 2 channels in input image.')

metadata = imgdata{1, 2};
LUT = imgdata{1, 3};
OME = imgdata{1, 4};

px_dist = str2double(metadata.get('Global CameraFramePixelDistance 0'));

ch1_pixel_data = imgdata{1, 1}{1, 1};
ch1_img = double(ch1_pixel_data)/double(max(ch1_pixel_data(:)));

ch2_pixel_data = imgdata{1, 1}{2, 1};
ch2_img = double(ch2_pixel_data)/double(max(ch2_pixel_data(:)));

nrows = size(ch1_img, 1);
ncols = size(ch1_img, 2);

composite_img = zeros(nrows, ncols, 3);
composite_img(:,:,2) = ch2_img;
composite_img(:,:,3) = ch1_img;

figure;
enter_pressed = false;
imshow(composite_img)
axis equal
hold on
i = 1;
while(enter_pressed == false)
    [x, y] = ginputWhite(1);
    try
        coordinates(i,:) = [x, y];
    catch ME
        break
    end
    if i == 1
        sData = scatter(coordinates(:,1), coordinates(:,2), 40, 'r');
    else
        delete(sData)
        plot(coordinates(:,1), coordinates(:,2), 'ro-', 'MarkerSize', 6);
    end
    i = i + 1;
end
spcv = cscvn([coordinates(:,1), coordinates(:,2)].');
hold off
close

Ltotal = arclength(coordinates(:,1), coordinates(:,2), 'spline');
distlist = (0:5:Ltotal)/Ltotal;
xyhat = interparc(distlist, coordinates(:,1), coordinates(:,2), 'spline');

width_chosen = false;
width = 20;

figure;

while width_chosen == false
    
    [x_inner, y_inner, x_outer, y_outer, R, unv, concavity, overlap] = parallel_curve(xyhat(:, 1), xyhat(:,2), width, 0, 0);
    
    imshow(composite_img)
    hold on
    axis equal
    fnplt(spcv, 'r--')
    plot(x_inner, y_inner, 'y--');
    plot(x_outer, y_outer, 'y--');
    hold off
    
    answer = questdlg('Would you like to change the width?', 'Band Width', 'Yes', 'No', 'No');
    switch answer
        case 'Yes'
            % prompt new width
            usr_input = inputdlg({'Enter new width:'}, 'Band Width', [1 35], {num2str(width)});
            width = str2double(char(usr_input));
            clf('reset')
        case 'No'
            width_chosen = true;
    end
    
end
img_file_parts = split(imgfile, '.');
output_name = strcat(char(img_file_parts{1,1}), '_ROI.fig');
savefig(output_name);

x_lb = ceil(min(min(x_inner), min(x_outer)));
x_ub = floor(max(max(x_inner), max(x_outer)));
y_lb = ceil(min(min(y_inner), min(y_outer)));
y_ub = floor(max(max(y_inner), max(y_outer)));

k_inner = dsearchn([x_inner, y_inner], xyhat);
k_outer = dsearchn([x_outer, y_outer], xyhat);

gray_img = zeros(nrows, ncols);
gray_img(y_lb:y_ub,x_lb:x_ub) = ch2_img(y_lb:y_ub,x_lb:x_ub);
figure
imshow(gray_img)
hold on
plot([x_lb, x_lb], [y_lb, y_ub], 'r-');
plot([x_lb, x_ub], [y_lb, y_lb], 'r-');
plot([x_ub, x_lb], [y_ub, y_ub], 'r-');
plot([x_ub, x_ub], [y_ub, y_lb], 'r-');

mean_intensities = zeros(length(xyhat(:,1)), 1);
var_intensities = zeros(length(xyhat(:,1)), 1);
for cnt = 1 : length(xyhat(:,1))
    xi = x_inner(k_inner(cnt));
    yi = y_inner(k_inner(cnt));
    xo = x_outer(k_outer(cnt));
    yo = y_outer(k_outer(cnt));
    ln = plot([xi, xo], [yi, yo], 'g-');
    ln.Color(4) = 0.25;
    q = improfile(ch2_pixel_data, [xi, xo], [yi, yo]);
    mean_intensities(cnt) = mean(q);
    var_intensities(cnt) = std(q);
end    
xlim([1, ncols]);
ylim([1, nrows]);
output_name = strcat(char(img_file_parts{1,1}), '_BBox.fig');
savefig(output_name);

figure
xvals = distlist.*Ltotal;
plot(xvals, mean_intensities, 'k-', 'LineWidth', 1.5)
hold on
opts = {'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.3, 'EdgeAlpha', 0.0};
fill_between(xvals, mean_intensities-0.5.*var_intensities, mean_intensities+0.5.*var_intensities, [], opts{:});
xlim([min(xvals)-0.05*(max(xvals)-min(xvals)), max(xvals)+0.05*(max(xvals)-min(xvals))])
xlabel('Length Along Curve (pixels)')
ylabel('Image Intensity')
output_name = strcat(char(img_file_parts{1,1}), '_Intensity.fig');
savefig(output_name);