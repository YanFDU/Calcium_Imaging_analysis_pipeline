clear all;
load('Fall.mat');

frame_rate = 30;
min_prominence = 0.4;
smooth_window = 5;
min_distance_sec = 1;
amplitude_scale_factor = 10;

accepted_cells = find(iscell(:,1) == 1);
mean_decay_times = [];

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
    peak_locs = peak_locs(valid);
    peak_vals = peak_vals(valid);

    
    min_dist_frames = round(min_distance_sec * frame_rate);
    final_locs = [];
    final_vals = [];
    last = -Inf;
    for i = 1:length(peak_locs)
        if isempty(final_locs) || (peak_locs(i) - last) >= min_dist_frames
            final_locs(end+1) = peak_locs(i);
            final_vals(end+1) = peak_vals(i);
            last = peak_locs(i);
        end
    end

    
    decay_t = [];
    for i = 1:length(final_locs)
        peak_idx = final_locs(i);
        peak_val = final_vals(i);
        decay_threshold = peak_val * exp(-1);

        for j = peak_idx+1:length(dff)
            if dff(j) <= decay_threshold
                decay_t(end+1) = (j - peak_idx) / frame_rate;
                break;
            end
        end
    end

    if ~isempty(decay_t)
        mean_decay_times(end+1) = mean(decay_t);
    end
end


figure;
histogram(mean_decay_times, 20, 'FaceColor', [1 0.6 0.2]);
hold on;
[f, xi] = ksdensity(mean_decay_times);
plot(xi, f * max(histcounts(mean_decay_times,20))/max(f), 'r-', 'LineWidth', 2);

m_decay = mean(mean_decay_times);
xline(m_decay, '--k', sprintf('Mean: %.2f s', m_decay), 'LabelOrientation','horizontal');

xlabel('Mean Decay Time (s)');
ylabel('Cell Count');
title('Decay Time Distribution');
grid on;
