
close all

% Define parameters
ppm_axis = linspace(cppm - ppm/2, cppm + ppm/2, size(rc123,1)); % Frequency axis in ppm
peak_ppm = [160.28, 163, 170.68, 183.16]; % Peaks of interest
%peak_ppm = [160.28, 163, 170.68]; % Peaks of interest
%peak_ppm = [148, 166, 184, 197]; % Peaks of interest
num_frames = 15; % Number of time frames

% Find nearest indices for each peak
[~, peak_idx] = min(abs(ppm_axis' - peak_ppm), [], 1);


y_labels = ["Peak 1", "Peak 2", "Peak 3", "Peak 4"];
y_positions = round(linspace(1, 8, 4));

% Normalize igray (proton anatomical image)
igray = double(igray);
igray = (igray - min(igray(:))) / (max(igray(:)) - min(igray(:))); % Normalize

% Preallocate an image array for montage
montage_images = zeros(size(igray,1), size(igray,2), 3, num_frames * length(peak_ppm));

% Create overlay images for each time frame and ppm peak
frame_idx = 1;
grid_thickness = 0; % Set grid line thickness (increase for thicker lines)


for p = 1:length(peak_ppm)

    for t = 1:num_frames
%         % Extract 8x8 CSI data at each peak
%         heatmap_data = abs(squeeze(rc123(peak_idx(p), :, :, t))); 
%         
%         % Resize CSI data to match igray dimensions
%         csi_resized = imresize(heatmap_data, [size(igray,1), size(igray,2)], 'lanczos3');
%         
%         % Apply colormap to CSI image
%         %csi_mapped = ind2rgb(gray2ind(mat2gray(csi_resized), 256), jet(256));
%         csi_mapped = ind2rgb(gray2ind(mat2gray(csi_resized, [0 50]), 256), jet(256));
% 
%         
%         % Overlay CSI onto igray image
%         alpha = 0.3; % Transparency factor
%         overlay_img = (1-alpha) * repmat(igray, [1, 1, 3]) + alpha * csi_mapped;
%         
%         % Store the overlay image for montage
%         montage_images(:,:,:,frame_idx) = overlay_img;
%         frame_idx = frame_idx + 1;
% 
%         % Add Y-axis labels at specific positions
%         if any(y == y_positions) && x == 1 % Apply only to the first column
%             ylabel(y_labels(y_positions == y), 'y','FontSize', 10, 'FontWeight', 'bold');
%         end
 
        % Extract 8x8 CSI data at each peak
        heatmap_data = abs(squeeze(rc123(peak_idx(p), :, :, t))); 
        
        % Resize CSI data (20% smaller than igray)
        new_size = round(0.64 * size(igray)); % Scale down by 20%
        csi_resized = imresize(heatmap_data, new_size, 'lanczos3');

        % Define fixed intensity range for color scaling
        csi_min = 0;  % Minimum expected value
        csi_max = 50; % Maximum expected value

        % Clip and normalize CSI data for consistent colormap mapping
        csi_resized_clipped = max(min(csi_resized, csi_max), csi_min);
        csi_resized_norm = (csi_resized_clipped - csi_min) / (csi_max - csi_min);

        % Apply colormap to CSI image
        csi_mapped = flipud(rot90(fliplr(ind2rgb(gray2ind(csi_resized_norm, 256), jet(256)))));

        % Ensure igray is in RGB format for overlay
        igray_rgb = repmat(igray, [1, 1, 3]);

        % Compute the overlay position (isocenter)
        igray_size = size(igray);
        csi_size = size(csi_mapped, 1:2); % Get 2D size of csi_mapped

        x_offset = round((igray_size(1) - csi_size(1)) / 2); % Center X
        y_offset = round((igray_size(2) - csi_size(2)) / 2); % Center Y

        % ------ Add Thin 8x8 Grid on the Heatmap ------
        grid_color = reshape([1, 1, 1], [1, 1, 3]); % Black grid lines

        % Define grid spacing based on new size
        grid_x = round(linspace(1, csi_size(2), 9)); % 9 points to create 8 columns
        grid_y = round(linspace(1, csi_size(1), 9)); % 9 points to create 8 rows

      % Draw vertical grid lines with thickness
        for gx = grid_x
            if gx > 1 && gx < csi_size(2) % Avoid out-of-bounds indexing
                csi_mapped(:, max(1, gx-grid_thickness):min(csi_size(2), gx+grid_thickness), :) = ...
                    repmat(grid_color, [csi_size(1), 2*grid_thickness+1, 1]);
            end
        end

        % Draw horizontal grid lines with thickness
        for gy = grid_y
            if gy > 1 && gy < csi_size(1) % Avoid out-of-bounds indexing
                csi_mapped(max(1, gy-grid_thickness):min(csi_size(1), gy+grid_thickness), :, :) = ...
                    repmat(grid_color, [2*grid_thickness+1, csi_size(2), 1]);
            end
        end
        % ---------------------------------------------

        % Create an empty overlay image (same size as igray)
        overlay = igray_rgb;

        % Insert csi_mapped into overlay at computed position
        overlay(x_offset:x_offset+csi_size(1)-1, y_offset:y_offset+csi_size(2)-1, :) = ...
            (1 - 0.3) * overlay(x_offset:x_offset+csi_size(1)-1, y_offset:y_offset+csi_size(2)-1, :) + ...
            0.3 * csi_mapped; % Apply 30% transparency

        % Store the overlay image for montage
        montage_images(:,:,:,frame_idx) = overlay;
        frame_idx = frame_idx + 1;

        % Add Y-axis labels at specific positions
        if any(y == y_positions) && x == 1 % Apply only to the first column
            ylabel(y_labels(y_positions == y), 'FontSize', 10, 'FontWeight', 'bold');
        end
    end
end

% Display montage (4 rows for peaks, 20 columns for time frames)
figure;
montage(montage_images, 'Size', [length(peak_ppm), num_frames], 'ThumbnailSize', [] );
title('Dynamic CSI Heatmaps (Peak nymber × Time Frames)');


%%
% Compute magnitude (since rc1 is complex)
rc1_mag = abs(rc1);

% Define ppm axis
ppm_axis = linspace(cppm - ppm/2, cppm + ppm/2, size(rc1, 1));

% Initialize storage for peak data
%top_peaks = struct('Intensity', [], 'PPM', [], 'X', [], 'Y', [], 'Time', []);

% Initialize storage for peak data
top_peaks = struct('AUC', [], 'PPM', [], 'X', [], 'Y', [], 'Time', []);

% Loop over all voxels (8x8) and time frames (15)
for x = 5
    for y = 5
        for t = 1
            % Extract spectrum for the voxel at time t
            spectrum = squeeze(rc1_mag(:, x, y, t));
            
            % Find peaks in the spectrum
            [peak_intensities, peak_locs] = findpeaks(spectrum, ppm_axis, 'SortStr', 'descend');
            
            % Keep only the top 4 peaks
            num_peaks = min(4, length(peak_intensities));
            
            for p = 1:num_peaks
                % Find indices within ±1 ppm range
                ppm_range = (ppm_axis >= peak_locs(p) - 1) & (ppm_axis <= peak_locs(p) + 1);
                
                % Compute the area under the curve (AUC) using trapezoidal integration
                auc_value = trapz(ppm_axis(ppm_range), spectrum(ppm_range));
                
                % Store peak details
                top_peaks(end+1) = struct(...
                    'AUC', auc_value, ...
                    'PPM', peak_locs(p), ...
                    'X', x, ...
                    'Y', y, ...
                    'Time', t ...
                );
            end
        end
    end
end


% Loop over all voxels (8x8) and time frames (15)
for x = 5
    for y = 5
        for t = 1
            % Extract spectrum for the voxel at time t
            spectrum = squeeze(rc1_mag(:, x, y, t));
            
            % Find peaks in the spectrum
            [peak_intensities, peak_locs] = findpeaks(spectrum, ppm_axis, 'SortStr', 'descend');
            
            % Keep only the top 4 peaks
            num_peaks = min(4, length(peak_intensities));
            
            for p = 1:num_peaks
                % Store peak details
                top_peaks(end+1) = struct(...
                    'Intensity', peak_intensities(p), ...
                    'PPM', peak_locs(p), ...
                    'X', x, ...
                    'Y', y, ...
                    'Time', t ...
                );
            end
        end
    end
end

% Remove the first empty entry
top_peaks(1) = [];

% Sort all detected peaks by intensity (descending)
[~, sort_idx] = sort([top_peaks.Intensity], 'descend');
top_peaks = top_peaks(sort_idx);

% Display the highest 4 peaks overall
fprintf('Top 4 Peaks in rc1:\n');
for i = 1:min(4, length(top_peaks))
    fprintf('Peak %d: Intensity = %.4f, PPM = %.2f, Voxel = (%d,%d), Time Frame = %d\n', ...
        i, top_peaks(i).Intensity, top_peaks(i).PPM, top_peaks(i).X, top_peaks(i).Y, top_peaks(i).Time);
end