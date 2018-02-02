%
% Segment and quantify cell clusters
% Channel 1: nuclei (RFP), Channel 2: cytoplasm (GFP)
% Imaging time: 60 hrs (15 min intervals), 240 frames
% Image properties: 0.32 um/px, 12 bit, 2560 X 2160 pixels, 20 X zoom
% Authors: Dhananjay Bhaskar <dhananjay_bhaskar@brown.edu> and Dr. Ian Wong <ian_wong@brown.edu>
% Last Modified: Jan 31, 2018
%

function [] = ClusterSegment()
    
    close all;

    nucleiImgs = dir(strcat('egf_XY1C1T', '*.tif'));
    cytoplasmImgs = dir(strcat('egf_XY1C2T', '*.tif'));

    if (length(nucleiImgs) ~= length(cytoplasmImgs))
        disp('Number of images does not match.');
        disp(length(nucleiImgs));
        disp(length(cytoplasmImgs));
        return;
    end
    
    % Get well number
    cDir = strsplit(pwd, filesep);
    cDir = cDir{end};
    well = cDir(1:3);
    
    % Create CSV file to store results
    header = {'Frame', 'Well', 'Cluster ID', 'Nuclei' 'CentroidX', 'CentroidY', 'Area', 'Perimeter', ...
             'Extent', 'Solidity', 'Compactness', 'Convexity', 'Circularity', 'Equivalent Circle Diameter', ...
             'Major Axis Length', 'Minor Axis Length', 'Eccentricity', 'Radius of gyration', ...
             'Aspect Ratio', 'Fractal Dimension', 'Elongation', 'Fractional Anisotropy'};
    
    CSVfile = strcat('MATLAB_Segmentation_', well ,'.csv');
    fid = fopen(CSVfile, 'w');
    fprintf(fid, '%s,', header{1,1:end-1});
    fprintf(fid, '%s\n', header{1,end});
    fclose(fid);
    
    % Create directories to store results
    mkdir('BW_Mask');
    mkdir(strcat('BW_Mask', filesep, 'Nuclei'));
    mkdir(strcat('BW_Mask', filesep, 'Cytoplasm'));
    mkdir('Watershed_Overlap');
    mkdir('Segmented_Boundary');
    mkdir('Betti_Intervals');
    mkdir('Persistence_Homology');
    
    for i = 1 : length(nucleiImgs)

        % Load images
        [~, Inucleifname, ~] = fileparts(nucleiImgs(i).name);
        Inucleitime = str2double(Inucleifname(end-2:end));
        Inuclei = imread(nucleiImgs(i).name);

        [~, Icytoplasmfname, ~] = fileparts(cytoplasmImgs(i).name);
        Icytoplasmtime = str2double(Icytoplasmfname(end-2:end));
        Icytoplasm = imread(cytoplasmImgs(i).name);
        
        % Get time
        if (Inucleitime ~= Icytoplasmtime)
            warning("Error: Images numbers do not match!");
            exit();
        end
        t_string = sprintf('%03d', Inucleitime);
        
        % Create folder for cell images
        mkdir(strcat('egf_', well, '_T', t_string));
        mkdir(strcat('Graph_', well, '_T', t_string));
        
        % Convert to grayscale
        if size(Inuclei, 3) == 3
            Inuclei = rgb2gray(Inuclei);
        end

        if size(Icytoplasm, 3) == 3
            Icytoplasm = rgb2gray(Icytoplasm);
        end
        
        % Display flag
        dp = 1;
        
        [Inuclei_processed, BWnuclei, BW_nuclei_overlap] = preprocess_nuclei(Inuclei, dp);
        imwrite(BW_nuclei_overlap, strcat('BW_Mask', filesep, 'Nuclei', filesep, 'egf_', well, '_T', t_string, '.png'), 'png');
        
        [Icytoplasm_processed, BWcytoplasm, BW_overlap] = preprocess_cytoplasm(Icytoplasm, BWnuclei, dp);
        imwrite(BW_overlap, strcat('BW_Mask', filesep, 'Cytoplasm', filesep, 'egf_', well, '_T', t_string, '.png'), 'png');
        
        [segmented_cells, segmented_borders, w_overlap] = watershed_segment(Icytoplasm_processed, BWcytoplasm, dp);
        imwrite(w_overlap, strcat('Watershed_Overlap', filesep, 'egf_', well, '_T', t_string, '.png'), 'png');
        
        % Extract features and write pif file
        feats = extract_features(t_string, well, segmented_cells, segmented_borders, BWnuclei, dp);
        dlmwrite(CSVfile, feats, '-append', 'precision', '%.3f', 'delimiter', ',');
        
        % pause before processing next image
        if dp == 1
            pause
        end
        
    end
    
end


function [res] = grayscale_remap(image, percentile)

    prc = prctile(image(:), percentile);
    prc = double(prc);
    map = gray(prc);
    res = ind2gray(image, map);
    res = im2double(res);
    
end


function [res, bw, bw_overlap] = preprocess_nuclei(image, dp)

    res = medfilt2(grayscale_remap(medfilt2(image), 99));
    
    bw = imbinarize(res, 'adaptive', 'Sensitivity', 0.4);

    % clear nuclei smaller than 21 um^2
    bw = bwareaopen(bw, 205);
    
    bw = imfill(bw, 'holes');
    
    bw = imopen(bw, strel('disk', 5));
    
    bw_overlap = res;
    bw_overlap(imdilate(bwperim(bw), ones(3, 3))) = 255;
    
    if (usejava('desktop') == 1 && dp == 1)
        figure
        subplot(1,2,1), imshow(res), title('Nuclei (RFP)')
        subplot(1,2,2), imshow(bw), title('Binary Foreground')
    end
    
end


function [res, bw, bw_overlap] = preprocess_cytoplasm(image, BWnuclei, dp)

    I = grayscale_remap(medfilt2(image), 99.9);

    background = imopen(I, strel('disk', 200));
    
    h = fspecial('gaussian', 10, 0.9);
    background = imfilter(background, h, 'replicate');
    
    nuclei_seg = regionprops(BWnuclei, 'Centroid');
    nuclei_centroids = cat(1, nuclei_seg.Centroid);
    nuclei_centroids = round(nuclei_centroids);

    res = I - background;
    res_bin = imadjust(res);
    res(res < 0.0) = 0;
    res_bin(res_bin < 0.0) = 0;
    
    BWnuclei_closed = imclose(BWnuclei, strel('disk', 100));  
    BWnuclei_closed = imclearborder(BWnuclei_closed);
    bw = activecontour(res_bin, BWnuclei_closed, 1000, 'Chan-Vese');
    
    bw = medfilt2(bw);
    
    bw = imclose(bw, strel('disk', 50));
    bw = imfill(bw, 'holes');
    
    % clear cells/clusters smaller than 256 um^2
    bw = bwareaopen(bw, 2500);
    
    L = bwlabel(bw, 4);
    [rows, cols] = size(L);
    bw = zeros(rows, cols);
    
    for k = 1 : size(nuclei_centroids)
        if isnan(nuclei_centroids(k,1))
            continue
        end
        x = nuclei_centroids(k,1);
        y = nuclei_centroids(k,2);
        if L(y, x) > 0
            val = L(y, x);
            bw(L == val) = 1;
        end
    end
    
    bw_overlap = res;
    bw_overlap(imdilate(bwperim(bw), ones(3, 3))) = 255;

    if (usejava('desktop') == 1 && dp == 1)
        figure
        subplot(2,2,1), imshow(I), title('Cytoplasm (GFP)')
        subplot(2,2,2), imshow(background), title('Estimated Background')
        subplot(2,2,3), imshow(res), title('Image - Background')
        subplot(2,2,4), imshow(bw), title('Binary Foreground')
    end

end


function [w_cells, w_borders, w_overlap] = watershed_segment(image, bin, dp)
    
    sigma = 2;
    smoothimage = imgaussfilt(image, sigma);
    [gradmag, ~] = imgradient(smoothimage, 'CentralDifference');
    
    % background markers
    D = bwdist(bin);
    DL = watershed(D);
    bgm = DL == 0;
    
	watershed_ridge = image;
	watershed_ridge(imdilate(bgm, ones(3, 3))) = 255;
    
    gradmag2 = imimposemin(gradmag, bgm | bin);
    L = watershed(gradmag2);
    
    w_overlap = image;
    
    classes = numel(unique(reshape(L, [1, numel(L)])));
	w_borders = zeros(size(L, 1), size(L, 2));
	w_cells = zeros(size(L, 1), size(L, 2));
    
    for i = 1 : classes
		tmp = zeros(size(L, 1), size(L, 2));
		tmp(L == i) = 1;
		w_cells(L == i) = i;
		w_borders(bwperim(tmp) > 0) = i;
    end

    % Remove segmentations overlapping border
    w_cells(ismember(w_cells, union(w_cells(:, [1 end]), w_cells([1 end], :)))) = 0;
    w_borders(ismember(w_borders, union(w_borders(:, [1 end]), w_borders([1 end], :)))) = 0;
    
    w_overlap(imdilate(w_borders > 0, ones(3, 3))) = 255;
    
    if (usejava('desktop') == 1 && dp == 1)
		figure
		subplot(2,2,1), imagesc(gradmag), colorbar, title('Gradient Magnitude')
        subplot(2,2,2), imshow(bin), title('Foreground Markers')
		subplot(2,2,3), imshow(watershed_ridge), title('Background Markers')
		subplot(2,2,4), imshow(w_overlap), title('Segmented Cell Boundary')
    end
    
end


function [result] = extract_features(t_string, well, segmented_cells, segmented_borders, BWnuclei, dp)

    % Add external libraries to path
    javaaddpath('./lib/javaplex.jar');
    import edu.stanford.math.plex4.*;
    javaaddpath('./lib/plex-viewer.jar');
    import edu.stanford.math.plex_viewer.*;

    cd './utility';
    addpath(pwd);
    cd '..';

    % Compute features from cytoplasm segmentation
    seg = regionprops(segmented_cells, 'Centroid', 'BoundingBox', 'Area', 'Perimeter', 'Extent', 'Solidity', ...
    'EquivDiameter', 'MajorAxisLength', 'MinorAxisLength', 'ConvexImage', 'Eccentricity');
    
    centroids = cat(1, seg.Centroid);
    bboxes = cat(1, seg.BoundingBox);
    area = [seg.Area];
    perim = [seg.Perimeter];
    extent = [seg.Extent];
    solidity = [seg.Solidity];
    eqDiameter = [seg.EquivDiameter];
    maj_axis = [seg.MajorAxisLength];
    min_axis = [seg.MinorAxisLength];
    ecc = [seg.Eccentricity];
    
    % Compute convex hull perimeter and fractal dimension for each cell
    convex_perimeter = zeros(size(seg, 1), 1);
    fractal_dimension = zeros(size(seg, 1), 1);
    
    for j = 1 : size(seg)
        
        tmp = zeros(size(segmented_cells, 1), size(segmented_cells, 2));
        tmp(segmented_cells == j) = 1;
        
        struct_array = regionprops(tmp, 'centroid', 'conveximage', 'boundingbox');
        
        if ~isempty(struct_array) && isfield(struct_array, 'Centroid')
            if ~isempty(struct_array.Centroid) && numel(extractfield(struct_array, 'Centroid')) == 2
                
                convex_hull = struct_array.ConvexImage;
                rect = struct_array.BoundingBox;
                struct_array = regionprops(convex_hull, 'perimeter');
                convex_perimeter(j) = struct_array.Perimeter;
                
                xmin = max(0, floor(rect(2)-1));
                xmax = min(size(segmented_cells, 1), floor(rect(2)+rect(4)+1));
                ymin = max(0, floor(rect(1)-1));
                ymax = min(size(segmented_cells, 2), floor(rect(1)+rect(3)+1));
                
                subimg = tmp(xmin:xmax,ymin:ymax);
                [n, r] = boxcount(subimg);
                p = polyfit(log(r), log(n), 1);
                beta_1 = p(1);
                lm = r.^beta_1 .* exp(p(2));
                fractal_dimension(j) = -beta_1;
                D = -gradient(log(n))./gradient(log(r));
                
                if (usejava('desktop') == 1 && dp == 1)
                    xls = 'r, Box Size';
                    yls = 'n(r), # Boxes';
                    figure
                    subplot(1,3,1), imshow(subimg)
                    subplot(1,3,2), semilogx(r, D, 'bo-'), xlabel(xls), ylabel('D, Local Fractal Dimension')
                    subplot(1,3,3)
                        loglog(r, n, 'bo-', r, (r/r(end)).^(-2), 'r--'), hold on
                        loglog(r, lm, 'g--'), xlabel(xls), ylabel(yls), hold off
                    legend('Box Count', 'Space-filling Count', 'Regression', 'Location', 'best');
                    title(strcat('OLS \beta_1 = ', string(beta_1)))
                    
                    sfname = strcat('egf_', well, '_T', t_string, filesep, 'Cluster_', int2str(j), '.png');
                    print(sfname, '-dpng');
                end
                
            end
        end
        
    end
    
    % Calculate features
    rgy = sqrt(maj_axis.^2 + min_axis.^2);
    aspect_ratio = maj_axis./min_axis;
    elongation = (maj_axis.^2 - min_axis.^2)./(maj_axis.^2 + min_axis.^2);
    frac_anisotropy = (maj_axis.^2 - min_axis.^2)./sqrt(maj_axis.^4 + min_axis.^4);
    compactness = sqrt((4*area)/pi)./maj_axis;
    convexity = convex_perimeter./perim;
    circularity = (4*pi*area)./(perim.^2);
        
    nuclei_seg = regionprops(BWnuclei, 'Centroid');
    nuclei_centroids = cat(1, nuclei_seg.Centroid);
    
    % Perform TDA (code borrowed from Dr. Ian Wong)
    max_dimension = 2;
    max_filtration_value = 200;
    num_divisions = 1000;
    
    % Find clusters using nuclei positions and plot connectivity graph
    min_thres_val = 100;
    max_thres_val = 200;
    cluster_list = find_clusters(nuclei_centroids, min_thres_val, max_thres_val, well, t_string);
    
    % Generate Vietoris-Rips stream
    stream = api.Plex4.createVietorisRipsStream(nuclei_centroids, max_dimension, max_filtration_value, num_divisions);
    % Calculate persistence
    persistence = api.Plex4.getModularSimplicialAlgorithm(max_dimension, 2);
    % Calculate intervals
    intervals = persistence.computeIntervals(stream);
    
    % Plot Betti intervals
    options.filename = strcat('Betti-egf-', well, '-T', t_string, '.png');
    options.max_filtration_value = max_filtration_value;
    options.max_dimension = max_dimension - 1;
    figH = plot_barcodes(intervals, options);
    set(figH, 'Name', strcat('Barcodes ', t_string, ' frame'), 'NumberTitle', 'off')
    system('mv Betti-*.png Betti_Intervals/')
    
    % Generate persistence plot
    figure('Name', strcat('Persistence ', t_string,' frame'), 'NumberTitle', 'off')
    % Betti 0
    startendpts0 = homology.barcodes.BarcodeUtility.getEndpoints(intervals, 0, false);
    % Betti 1
    startendpts1 = homology.barcodes.BarcodeUtility.getEndpoints(intervals, 1, false);
    plot(startendpts0(:,1), startendpts0(:,2), 'ro'); hold on
    % If 2-simplex is present                                           
    if ~isempty(startendpts1)
        plot(startendpts1(:,1), startendpts1(:,2), 'bo')
    end
    plot(0:200, 0:200, 'k-')
    xlabel('Start \epsilon')
    ylabel('End \epsilon') 
    hold off
    print(strcat('Persistence_Homology', filesep, 'egf_', well, '_T', t_string, '.png'), '-dpng');
                                            
    w_borders_color = label2rgb(imdilate(segmented_borders, ones(3, 3)), 'jet', [.9 .9 .9], 'shuffle');
    
    for k = 1 : size(nuclei_seg)
        if isnan(nuclei_centroids(k,1))
            continue
        end
        w_borders_color = insertMarker(w_borders_color, [nuclei_centroids(k,1) nuclei_centroids(k,2)], '*', 'color', 'blue', 'size', 15);
    end
    
    % Assign cluster/cell ID, display segmented boundary and count nuclei
    if (usejava('desktop') == 1 && dp == 1)

        figure
        imshow(w_borders_color), title('Watershed (Marker) Outlines')
        hold on

        for j = 1 : size(seg)
            
            if isnan(centroids(j,1))
                continue
            end
            
            tmp = zeros(size(segmented_cells, 1), size(segmented_cells, 2));
            tmp(segmented_cells == j) = 1;
            
            struct_array = regionprops(tmp, 'centroid', 'boundingbox');
            if ~isempty(struct_array) && isfield(struct_array, 'Centroid')
                if ~isempty(struct_array.Centroid) && numel(extractfield(struct_array, 'Centroid')) == 2
                    
                    centroid = cat(1, struct_array.Centroid);
                    rect = struct_array.BoundingBox;
                    num_nuclei = 0;
                    
                    for k = 1 : size(nuclei_seg)
                        if isnan(nuclei_centroids(k,1))
                            continue
                        end
                        nx = nuclei_centroids(k,1);
                        ny = nuclei_centroids(k,2);
                        if nx > rect(1) && nx < rect(1) + rect(3)
                            if ny > rect(2) && ny < rect(2) + rect(4)
                                num_nuclei = num_nuclei + 1;
                            end
                        end
                    end
                    
                    display_string = strcat('ID: ', int2str(j), ', Nuclei: ', int2str(num_nuclei));
                    txt = text(centroid(1), centroid(2), display_string);
                    set(txt, 'fontsize', 10);
                    
                end
            end
            
        end
        print(strcat('Segmented_Boundary', filesep, 'egf_', well, '_T', t_string, '.png'), '-dpng');

    end
    
    % Count number of nuclei inside each cell (again) and collate features
    cnt = 1;
    result = [];
    segmentations = zeros(size(segmented_cells,1), size(segmented_cells,2));
    
    for j = 1 : size(seg)
        if isnan(centroids(j,1))
            continue
        end
        num_nuclei = 0;
        for k = 1 : size(nuclei_seg)
            if isnan(nuclei_centroids(k,1))
                continue
            end
            nx = nuclei_centroids(k,1);
            ny = nuclei_centroids(k,2);
            if nx > bboxes(j,1) && nx < bboxes(j,1) + bboxes(j,3)
                if ny > bboxes(j,2) && ny < bboxes(j,2) + bboxes(j,4)
                    num_nuclei = num_nuclei + 1;
                end
            end
        end
        fnum = str2double(t_string);
        wnum = str2double(well(2:end));
        result(cnt,:) = [fnum wnum j num_nuclei centroids(j,1) centroids(j,2) area(j) perim(j) ...
                      extent(j) solidity(j) compactness(j), convexity(j), circularity(j), eqDiameter(j) ...
                      maj_axis(j) min_axis(j) ecc(j) rgy(j) aspect_ratio(j) fractal_dimension(j) ...
                      elongation(j) frac_anisotropy(j)];
        
        segmentations(segmented_cells == j) = j;         
        cnt = cnt + 1;
    end
    
    % Serialize segmented cells in pif file for further processing
    write_pif_file(segmentations, t_string, well);
    num_clusters = cnt - 1;
    disp(strcat('Number of clusters saved: ', int2str(num_clusters)));
    
end


function [] = write_pif_file(segmentation_mat, t_string, well)

    cluster_ids = unique(nonzeros(segmentation_mat));
	numPIFrows = nnz(segmentation_mat);

	PIF_data = zeros(numPIFrows, 5);
	ind = 1;
    for i = 1 : length(cluster_ids)
		[rows, cols] = find(segmentation_mat == cluster_ids(i));
        for cnt = 1 : length(rows)
			PIF_data(ind, 1) = cluster_ids(i);
			PIF_data(ind, 2) = rows(cnt);
			PIF_data(ind, 3) = rows(cnt);
			PIF_data(ind, 4) = cols(cnt);
			PIF_data(ind, 5) = cols(cnt);
			ind = ind + 1;
        end
    end
    
    PIFname = strcat('egf_', well, '_T', t_string, filesep, 'Cell_Imgs.pif');
    fileID = fopen(PIFname, 'w');
    cnt = 1;
    
    while cnt < ind
        fprintf(fileID, '%d CellU %d %d %d %d\n', PIF_data(cnt,:));
        cnt = cnt + 1;
    end
    
    fclose(fileID);

end


% Attribution: Dr. Ian Wong
function [res] = find_clusters(XYCoords, minval, maxval, well, t_string)

    NumCol = maxval - minval + 30;

    % Initialize cluster
    storeClusters = repmat(1:size(XYCoords,1), 1, 1);
    storeClusters = storeClusters';

    % Store NaN locations & copy NaN locations to storeClusters
    storeNaN = ~isnan(XYCoords);
    storeClusters = storeClusters.*storeNaN;
    storeClusters(storeClusters==0) = NaN;

    storeAllNeighbors = nan(1,3);

    % Decrement search radius
    figure
    for threshold = minval:5:maxval

        counter = 1;

        % Calculate and generate storeClusters matrix
        storeNeighbors = nan(size(XYCoords,1), 3);

        % Compute distances and cluster cells if they are within threshold
        for i = 1:(size(XYCoords,1)-1)
            
            if isnan(XYCoords(i,1)) == 1
                continue
            end
            
            for j = (i+1):size(XYCoords,1)
                
                if isnan(XYCoords(j,1)) == 1
                    continue
                end

                x1 = XYCoords(i,1);
                x2 = XYCoords(j,1);
                y1 = XYCoords(i,2);
                y2 = XYCoords(j,2);

                Distance = sqrt((x2-x1)^2+(y2-y1)^2);

                if Distance <= threshold

                    storeNeighbors(counter,1) = i;
                    storeNeighbors(counter,2) = j;
                    storeNeighbors(counter,3) = threshold;

                    % what clusters do the two cells originally belong to?
                    Cluster_i = storeClusters(i,1);
                    Cluster_j = storeClusters(j,1);

                    % what other cells belong to Cluster_i and Cluster_j?
                    findCellsWithin_i = storeClusters(:,1)==Cluster_i;
                    findCellsWithin_j = storeClusters(:,1)==Cluster_j;

                    % rename all clusters into one ID
                    storeClusters(findCellsWithin_i,1) = Cluster_i;
                    storeClusters(findCellsWithin_j,1) = Cluster_i;
                    
                    counter = counter + 1;
                    
                end

            end
            
        end

        % Rename clusters in sequential order
        count = 1;
        for k = 1:size(storeClusters,1)
            
            findCluster = find(storeClusters(:,1) == k);
            
            if isempty(findCluster) == 1
                continue
            end
            storeClusters(findCluster,1) = count;
            count = count + 1;
            
        end

        % Find cells that are present at this time frame (not nan)
        findPresent = ~isnan(XYCoords(:,1));
        
        % Plot nuclei positions in red
        plot(XYCoords(findPresent,1), XYCoords(findPresent,2), 'ro', 'MarkerSize', 10);
        hold on
        
        % Plot cluster nuclei in blue
        for n = 1:count
            
            fTemp = find(storeClusters(:,1) == n);

            if length(fTemp) > 1
                plot(XYCoords(fTemp,1), XYCoords(fTemp,2), 'bo', 'MarkerSize', 10); 
                hold on
            end
            
        end

        flinks = ~isnan(storeNeighbors(:,1));
        storeNeighbors = storeNeighbors(flinks,:);

        cmapk = flipud(cbrewer('seq', 'Greys', NumCol, 'PCHIP'));
        
        % Plot edges between nuclei in clusters
        for l = 1:length(storeNeighbors)
            plot(XYCoords(storeNeighbors(l,1:2),1), XYCoords(storeNeighbors(l,1:2),2), ...
                '-', 'Color', cmapk(threshold-minval+1,:));
        end
        hold off
        axis equal;
        title(strcat('Number of cells/clusters found:', int2str(count)));
        set(gca, 'Ydir', 'reverse');
        drawnow;
        
        sfname = strcat('Graph_', well, '_T', t_string, filesep, 'Graph_Threshold_', int2str(threshold), '.png');
        print(sfname, '-dpng');
        
        storeAllNeighbors = [storeAllNeighbors; storeNeighbors];
        
        pause(1);

    end
    
    res = storeAllNeighbors(2:length(storeAllNeighbors),:);

end
