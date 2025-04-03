% Define parameters
ppm_axis = linspace(cppm - ppm/2, cppm + ppm/2, size(rc123,1)); % Frequency axis in ppm
peak_ppm = [150, 168, 189, 199]; % Peaks of interest
num_frames = 15; % Number of time frames

% Find nearest indices for each peak
[~, peak_idx] = min(abs(ppm_axis' - peak_ppm), [], 1);

% Normalize igray (proton anatomical image)
igray = double(igray);
igray = (igray - min(igray(:))) / (max(igray(:)) - min(igray(:))); % Normalize

% Preallocate an image array for montage
montage_images = zeros(size(igray,1), size(igray,2), 3, num_frames * length(peak_ppm));

% Create overlay images for each time frame and ppm peak
frame_idx = 1;
for p = 1:length(peak_ppm)
    for t = 1:num_frames
        % Extract 8x8 CSI data at each peak
        heatmap_data = abs(squeeze(rc123(peak_idx(p), :, :, t))); 
        
        % Resize CSI data to match igray dimensions with higher-quality interpolation
        csi_resized = imresize(heatmap_data, [size(igray,1), size(igray,2)], 'lanczos3');
        
        % Apply colormap with higher resolution (1024 levels instead of 256)
        csi_mapped = ind2rgb(gray2ind(mat2gray(csi_resized, [0 70]), 1024), jet(1024));

        % Overlay CSI onto igray image using high-quality fusion
        overlay_img = imfuse(repmat(igray, [1, 1, 3]), csi_mapped, 'blend', 'Scaling', 'joint');
        
        % Store the overlay image for montage
        montage_images(:,:,:,frame_idx) = overlay_img;
        frame_idx = frame_idx + 1;
    end
end

figure;
montage(montage_images, 'Size', [length(peak_ppm), num_frames],'ThumbnailSize',[]);
title('Dynamic CSI Heatmaps (4 Peaks × 15 Time Frames)');

% Display montage (4 rows for peaks, 15 columns for time frames) with high resolution
%figure;
%montage(montage_images, 'Size', [length(peak_ppm), num_frames], 'ThumbnailSize', []);
%title('Dynamic CSI Heatmaps (4 Peaks × 15 Time Frames)');

% Set high-resolution figure size for better visibility
set(gcf, 'Position', [100, 100, 2000, 1000]);

% Save the montage as a high-resolution image
%saveas(gcf, 'high_res_montage.png');