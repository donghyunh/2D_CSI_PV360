close all

% Define figure for montage
figure;
tiledlayout(8, 8, 'Padding', 'compact', 'TileSpacing', 'compact'); % 8x8 grid layout

% Define ppm axis (assuming you have 'cppm' and 'ppm' defined)
ppm_axis = linspace(cppm - ppm/2, cppm + ppm/2, size(rc1, 1)); % Frequency axis in ppm
max_intensity = max(abs(rc1(:)));

% Loop through each voxel (8x8 grid)
for x = 1:8
    for y = 1:8
        nexttile; % Create subplot in the grid
        
        % Extract spectrum (512 data points) across 15 time frames
        spectra = squeeze(rc1(:, x, y, 10)); % 512 × 15 matrix
        
        % Plot all 15 time frames with transparency
        plot(ppm_axis, flipud(abs(spectra)), 'LineWidth', 1.2); 
        %plot(ppm_axis, abs(spectra), 'LineWidth', 1.2, 'Color', [0 0 1 0.4]);
        hold on;
        
        % Overlay the mean spectrum in bold
        plot(ppm_axis, mean(flipud(abs(spectra)), 2), 'r', 'LineWidth', 1.5);
        
        
        % Formatting
        set(gca, 'XDir', 'reverse'); % Standard convention for ppm axis
        axis tight;
        ylim ([0 max_intensity/20]);
        %xticks([]); yticks([]); % Remove axis labels for a clean look
        title(sprintf('(%d,%d)', x, y), 'FontSize', 8); % Label voxel coordinates
    end
end

% Add a common title
sgtitle('Spectral Montage Across 8×8 Voxels', 'FontSize', 14);

%%
% Select a specific voxel
x = 6; % Change to the voxel of interest
y = 3;

% Extract spectrum (512 data points) across 15 time frames
spectra = squeeze(rc1(:, x, y, :))'; % 15 × 512 matrix (Time frames as rows)

% Define ppm axis (Ensure ppm_axis matches the data dimension)
%ppm_axis = linspace(10, -10, size(spectra, 2)); % Modify if needed

% Convert to magnitude spectrum
spectra = abs(spectra);

% Create a table with ppm values as column names
%T = array2table(spectra);

% Plot using stackedplot
figure;
h = stackedplot(spectra');

% Formatting
xlabel('Chemical Shift (ppm)');
ylabel('Signal Intensity');
axis off;
set(gcf, 'color', 'w')
title(sprintf('Stacked Spectra at Voxel (%d, %d)', x, y));
