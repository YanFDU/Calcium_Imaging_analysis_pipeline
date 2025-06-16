clear all;
load('Fall.mat');

frame_rate = 30;
amplitude_scale_factor = 10;
min_prominence = 0.4;
smooth_window = 5;
min_distance_sec = 1;

accepted_cells = find(iscell(:,1) == 1);
mean_amplitudes = [];

for idx = 1:length(accepted_cells)
    roi = accepted_cells(idx);
    F_trace = F(roi, :);
    Fneu_trace = Fneu(roi, :);
    r = 0.7;
    F_corr = F_trace - r * Fneu_trace;
    F_base = median(F_corr);
    dff = ((F_corr - F_base) ./ F_base) * amplitude_scale_factor;
    dff_smoothed = movmean(dff, smooth_window);

    
    is_peak = [false, dff_smoothed(2:end-1) > dff_smoothed(1:end-2) & dff_smoothed(2:end-1) > dff_smoothed(3:end), false];
    peak_locs = find(is_peak);
    peak_vals = dff_smoothed(peak_locs);
    valid = peak_vals > min_prominence;
    peak_vals = peak_vals(valid);
    peak_locs = peak_locs(valid);

    
    min_dist_frames = round(min_distance_sec * frame_rate);
    final_vals = [];
    last = -Inf;
    for i = 1:length(peak_locs)
        if isempty(final_vals) || (peak_locs(i) - last) >= min_dist_frames
            final_vals(end+1) = peak_vals(i);
            last = peak_locs(i);
        end
    end

    if ~isempty(final_vals)
        mean_amplitudes(end+1) = mean(final_vals);
    end
end


figure;
histogram(mean_amplitudes, 20, 'FaceColor', [0.2 0.6 1]);
hold on;


[f, xi] = ksdensity(mean_amplitudes);
plot(xi, f * max(histcounts(mean_amplitudes,20))/max(f), 'r-', 'LineWidth', 2);


m_amp = mean(mean_amplitudes);
xline(m_amp, '--k', sprintf('Mean: %.2f', m_amp), 'LabelOrientation','horizontal', 'LabelVerticalAlignment','bottom');

xlabel('Mean Peak Amplitude (ΔF/F scaled)');
ylabel('Cell Count');
title('Amplitude Distribution');
grid on;
