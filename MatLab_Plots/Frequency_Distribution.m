clear all;
load('Fall.mat');

frame_rate = 30;
amplitude_scale_factor = 10;  
min_prominence = 0.4;
smooth_window = 5;
min_distance_sec = 1;

accepted_cells = find(iscell(:,1) == 1);
peak_frequencies = [];

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

    
    min_dist_frames = round(min_distance_sec * frame_rate);
    cleaned = [];
    last = -Inf;
    for i = 1:length(peak_locs)
        if isempty(cleaned) || (peak_locs(i) - last) >= min_dist_frames
            cleaned(end+1) = peak_locs(i);
            last = peak_locs(i);
        end
    end

    
    if ~isempty(cleaned)
        duration_sec = length(dff) / frame_rate;
        freq = length(cleaned) / duration_sec;
        peak_frequencies(end+1) = freq;
    end
end


figure;
histogram(peak_frequencies, 20, 'FaceColor', [0.4 0.8 0.2]);
hold on;
[f, xi] = ksdensity(peak_frequencies);
plot(xi, f * max(histcounts(peak_frequencies,20))/max(f), 'r-', 'LineWidth', 2);

m_freq = mean(peak_frequencies);
xline(m_freq, '--k', sprintf('Mean: %.2f Hz', m_freq), 'LabelOrientation','horizontal');

xlabel('Peak Frequency (Hz)');
ylabel('Cell Count');
title('Frequency Distribution');
grid on;
