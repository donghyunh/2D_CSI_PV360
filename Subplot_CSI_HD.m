close all

% Define figure
figure;

% Loop through peaks (rows)
for p = 1:length(peak_ppm)

    % Determine min/max for each peak across time frames
    csi_min = 0;  % Fixed minimum
    csi_max = max(abs(rc123(peak_idx(p), :, :, :)), [], 'all'); % Peak-specific max

    % Loop through time frames (columns)
    for t = 1:num_frames
        subplot(length(peak_ppm), num_frames, (p-1)*num_frames + t);

        % Extract and resize heatmap data (scaled down to 64% of igray)
        heatmap_data = abs(squeeze(rc123(peak_idx(p), :, :, t))); 
        new_size = round(0.64 * size(igray)); % Reduce size by 20%
        csi_resized = imresize(heatmap_data, new_size, 'lanczos3');

        % Normalize and apply colormap
        csi_resized_norm = (max(min(csi_resized, csi_max), csi_min) - csi_min) / (csi_max - csi_min);
        csi_mapped = rot90(fliplr(ind2rgb(gray2ind(csi_resized_norm, 256), jet(256))));

        % Ensure igray is in RGB format
        igray_rgb = repmat(igray, [1, 1, 3]);

        % Compute offsets to center heatmap in igray
        igray_size = size(igray);
        csi_size = size(csi_mapped, 1:2); % Get 2D size of csi_mapped

        x_offset = round((igray_size(1) - csi_size(1)) / 2);
        y_offset = round((igray_size(2) - csi_size(2)) / 2);

        % Create an empty overlay image
        overlay = igray_rgb;

        % Insert the heatmap at computed position (centered)
        overlay(x_offset:x_offset+csi_size(1)-1, y_offset:y_offset+csi_size(2)-1, :) = ...
            (1 - 0.3) * overlay(x_offset:x_offset+csi_size(1)-1, y_offset:y_offset+csi_size(2)-1, :) + ...
            0.3 * csi_mapped; % Apply 30% transparency

        % Display the overlay image
        imagesc(overlay);
        axis off;
        colormap('jet');
        caxis([csi_min, csi_max]); % Individual scaling per row

        % Add title for the first row
        if p == 1
            title(['Frame ', num2str(t)]);
        end

        % Add labels for the first column
        if t == 1
            ylabel(['Peak ', num2str(p)], 'FontSize', 10, 'FontWeight', 'bold');
        end
    end

    % Add a colorbar for each row
    colorbar;
end

sgtitle('Dynamic CSI Heatmaps (Peaks × Time Frames)');
