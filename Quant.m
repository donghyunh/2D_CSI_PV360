close all;

% === Define Parameters ===
selected_voxels = [7,5; 7,4; 7,3]; % (Row, Column) positions
num_voxels = size(selected_voxels, 1);
num_frames = size(rc1, 4);
peak_ppm = [163.68, 170.68, 153.16]; % Peaks for fitting

% Define ppm axis
ppm_axis = linspace(cppm - ppm/2, cppm + ppm/2, size(rc1, 1));

% Find nearest indices for each peak
[~, peak_idx] = min(abs(ppm_axis' - peak_ppm), [], 1);

% Preallocate storage for AUC values
auc_curves = zeros(num_voxels, num_frames, length(peak_ppm));

% === Extract AUC for Each Voxel and Peak ===
for v = 1:num_voxels
    x = selected_voxels(v, 1);
    y = selected_voxels(v, 2);
    
    if x > size(rc1, 2) || y > size(rc1, 3)
        fprintf('Skipping voxel (%d, %d) - Out of bounds\n', x, y);
        continue;
    end

    for p = 1:length(peak_ppm)
        ppm_range = (ppm_axis >= peak_ppm(p) - 1) & (ppm_axis <= peak_ppm(p) + 1);
        
        if sum(ppm_range) < 2
            fprintf('Skipping peak at %.2f ppm - No valid range\n', peak_ppm(p));
            continue;
        end

        for t = 1:num_frames
            spectrum = abs(squeeze(rc1(:, x, y, t)));
            auc_curves(v, t, p) = trapz(ppm_axis(ppm_range), spectrum(ppm_range));
        end
    end
end

% === Plot AUC Dynamic Curves (Figure 1: Each subplot is a voxel with 3 peaks) ===
figure;
time_points = (1:num_frames)'; % Column vector
smooth_time = linspace(1, num_frames, 100)'; % More points for smoother curve
colors = lines(length(peak_ppm)); % Different colors for 160, 170, and 183 ppm peaks

for v = 1:num_voxels
    subplot(num_voxels, 1, v);
    hold on;

    for p = 1:length(peak_ppm)
        auc_values = squeeze(auc_curves(v, :, p))'; % Extract AUC values as row vector

        % Polynomial Regression Fit (Degree 2 for smooth curve)
        poly_coeffs = polyfit(time_points, auc_values, 4);
        smooth_fit = polyval(poly_coeffs, smooth_time);

        % Assign color for this peak
        peak_color = colors(p, :);

        % Plot Original Data Points
        scatter(time_points, auc_values, 60, 'filled', 'MarkerFaceColor', peak_color, ...
                'MarkerEdgeColor', 'k', 'DisplayName', sprintf('Peak %.2f ppm', peak_ppm(p)));

        % Plot Smooth Polynomial Fit
        plot(smooth_time, smooth_fit, '-', 'LineWidth', 2, 'Color', peak_color, ...
             'DisplayName', sprintf('PolyFit %.2f ppm', peak_ppm(p)));
    end

    hold off;
    %xlabel('Time Frame', 'FontSize', 12, 'FontWeight', 'bold');
    %ylabel('AUC Intensity', 'FontSize', 12, 'FontWeight', 'bold');
    %title(sprintf('Voxel (%d,%d): Dynamic Curves', selected_voxels(v,1), selected_voxels(v,2)), 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best');
    set(gcf, 'color', 'w')
    %axis off
    grid on;
end

% === Compute and Plot Dynamic Ratios (Figure 2: Each subplot is a voxel with 2 ratios) ===
figure;
ratio_labels = {'Bicarb/Pyr (160/170)', 'Lac/Pyr (183/170)'};
ratio_colors = [0 0.5 0; 0.8 0 0]; % Green for 160/170, Red for 183/170

for v = 1:num_voxels
    subplot(num_voxels, 1, v);
    hold on;

    for r = 1:2
        % Compute dynamic ratio while handling division by zero
        denom = squeeze(auc_curves(v, :, 2));
        denom(denom == 0) = NaN;  % Avoid division by zero
        
        if r == 1
            ratio_values = auc_curves(v, :, 1) ./ denom;
        else
            ratio_values = auc_curves(v, :, 3) ./ denom;
        end

        ratio_values = squeeze(ratio_values)'; % Convert to row vector

        % Polynomial Regression Fit for Ratios (Degree 2)
        poly_coeffs = polyfit(time_points, ratio_values, 2);
        smooth_fit = polyval(poly_coeffs, smooth_time);

        % Assign color for this ratio
        ratio_color = ratio_colors(r, :);

        % Plot Original Data Points
        scatter(time_points, ratio_values, 60, 'filled', 'MarkerFaceColor', ratio_color, ...
                'MarkerEdgeColor', 'k', 'DisplayName', ratio_labels{r});

        % Plot Smooth Polynomial Fit
        plot(smooth_time, smooth_fit, '-', 'LineWidth', 2, 'Color', ratio_color, ...
             'DisplayName', sprintf('PolyFit %s', ratio_labels{r}));
    end

    hold off;
    xlabel('Time Frame', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Dynamic Ratio', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Voxel (%d,%d): Ratio Dynamics', selected_voxels(v,1), selected_voxels(v,2)), 'FontSize', 14, 'FontWeight', 'bold');
    %legend('Location', 'best');
    set(gcf, 'color', 'w')
    grid on;
end
