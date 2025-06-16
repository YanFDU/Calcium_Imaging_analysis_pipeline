clear all;


load('Fall.mat');  


accepted_cells = iscell(:,1) == 1;
F = F(accepted_cells, :);
Fneu = Fneu(accepted_cells, :);
spks = spks(accepted_cells, :);


all_roi_indices = find(accepted_cells);


r = 0.7;
F_corrected = F - r * Fneu;


F_baseline = median(F_corrected, 2);
deltaF_F = (F_corrected - F_baseline) ./ F_baseline;


num_ROIs_to_plot = min(20, size(deltaF_F, 1));
roi_display_labels = all_roi_indices(1:num_ROIs_to_plot);  


num_frames = size(deltaF_F, 2);
time_in_seconds = (1:num_frames) / 30;  


amplitude_scale_factor = 300;
spks_scale_factor = 1;


deltaF_F = deltaF_F(1:num_ROIs_to_plot, :) * amplitude_scale_factor;
spks_to_plot = spks(1:num_ROIs_to_plot, :) * spks_scale_factor;


figure('Position', [100, 100, 400, 800]);
hold on;

offset = 0.3 * max(amplitude_scale_factor, spks_scale_factor);  % Reduced spacing
colors = lines(num_ROIs_to_plot);


for i = 1:num_ROIs_to_plot
    plot(time_in_seconds, spks_to_plot(i, :) + (i-1)*offset, 'Color', [0.3, 0.3, 0.3]);
end


for i = 1:num_ROIs_to_plot
    plot(time_in_seconds, deltaF_F(i, :) + (i-1)*offset, 'Color', colors(i,:) * 0.6, 'LineWidth', 1.5);
end


xlabel('Time (seconds)');
ylabel('ROI #');
title('ΔF/F & Deconvolved Traces');


set(gca, 'YTick', (0:offset:(num_ROIs_to_plot - 1)*offset));
set(gca, 'YTickLabel', roi_display_labels);


extra_margin = offset * 0.5;
ylim([-extra_margin, (num_ROIs_to_plot - 1)*offset + extra_margin]);

hold off;
