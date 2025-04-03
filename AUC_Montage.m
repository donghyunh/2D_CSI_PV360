close all;


% === Define Parameters ===
num_frames = 15; % Number of time frames
peak_ppm = [163.7, 171.86, 179.42, 183.3]; % Peaks of interest

% Define ppm axis based on center and range
ppm_axis = linspace(cppm - ppm/2, cppm + ppm/2, size(rc123, 1)); 

% Find nearest indices for each peak
[~, peak_idx] = min(abs(ppm_axis' - peak_ppm), [], 1);

% Normalize grayscale anatomical image
igray = double(igray);
igray = (igray - min(igray(:))) / (max(igray(:)) - min(igray(:))); 

% Preallocate montage storage
montage_images = zeros(size(igray,1), size(igray,2), 3, num_frames * length(peak_ppm));

% === Overlay CSI Heatmaps on Anatomical Image ===
frame_idx = 1;
grid_thickness = 0; % Set grid line thickness (increase for thicker lines)

for p = 1:length(peak_ppm)
    for t = 1:num_frames
        % Extract CSI data
        heatmap_data = abs(squeeze(rc123(peak_idx(p), :, :, t))); 
        
        % Resize CSI data
        new_size = round(0.8 * size(igray)); % Scale down by 20%
        csi_resized = imresize(heatmap_data, new_size, 'lanczos3');

        % Clip and normalize data
        csi_min = 0; csi_max = 300;
        csi_resized = max(min(csi_resized, csi_max), csi_min);
        csi_resized = (csi_resized - csi_min) / (csi_max - csi_min);

        % Apply colormap
        csi_mapped = flipud(rot90(fliplr(ind2rgb(gray2ind(csi_resized, 256), jet(256)))));

        % Compute overlay position
        igray_size = size(igray);
        csi_size = size(csi_mapped, 1:2);
        x_offset = round((igray_size(1) - csi_size(1)) / 2); 
        y_offset = round((igray_size(2) - csi_size(2)) / 2);

        % Create overlay
        igray_rgb = repmat(igray, [1, 1, 3]);
        overlay = igray_rgb;
        overlay(x_offset:x_offset+csi_size(1)-1, y_offset:y_offset+csi_size(2)-1, :) = ...
            (1 - 0.3) * overlay(x_offset:x_offset+csi_size(1)-1, y_offset:y_offset+csi_size(2)-1, :) + ...
            0.3 * csi_mapped;

        % Store for montage
        montage_images(:,:,:,frame_idx) = overlay;
        frame_idx = frame_idx + 1;
    end
end

% Display montage
figure;
montage(montage_images, 'Size', [length(peak_ppm), num_frames], 'ThumbnailSize', [] );
title('Dynamic CSI Heatmaps (Peak × Time Frames)');

% === Compute Area Under the Curve (AUC) for Each Peak ===
rc1_mag = abs(rc1); % Convert complex data to magnitude
ppm_axis = linspace(cppm - ppm/2, cppm + ppm/2, size(rc1, 1));

top_peaks = []; % Initialize storage

for x = 5
    for y = 5
        for t = 1
            spectrum = squeeze(rc1_mag(:, x, y, t)); % Extract spectrum
            
            % Find spectral peaks
            [peak_intensities, peak_locs] = findpeaks(spectrum, ppm_axis, 'SortStr', 'descend');
            num_peaks = min(4, length(peak_intensities));
            
            if isempty(peak_intensities), continue; end % Skip if no peaks found
            
            for p = 1:num_peaks
                % Define ±1 ppm range
                ppm_range = (ppm_axis >= peak_locs(p) - 1) & (ppm_axis <= peak_locs(p) + 1);
                
                if sum(ppm_range) < 2, continue; end % Ensure valid range
                
                % Compute AUC
                auc_value = trapz(ppm_axis(ppm_range), spectrum(ppm_range));
                
                % Store data
                top_peaks = [top_peaks; struct(...
                    'AUC', auc_value, ...
                    'PPM', peak_locs(p), ...
                    'X', x, ...
                    'Y', y, ...
                    'Time', t ...
                )];
            end
        end
    end
end

% Sort peaks by AUC
if isempty(top_peaks)
    fprintf('No peaks detected.\n');
else
    [~, sort_idx] = sort([top_peaks.AUC], 'descend');
    top_peaks = top_peaks(sort_idx);
    
    % Display top 4 peaks
    fprintf('Top 4 Peaks in rc1 (AUC-based):\n');
    for i = 1:min(4, length(top_peaks))
        fprintf('Peak %d: AUC = %.4f, PPM = %.2f, Voxel = (%d,%d), Time Frame = %d\n', ...
            i, top_peaks(i).AUC, top_peaks(i).PPM, top_peaks(i).X, top_peaks(i).Y, top_peaks(i).Time);
    end
end
