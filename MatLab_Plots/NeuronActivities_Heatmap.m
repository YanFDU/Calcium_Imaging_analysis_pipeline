clear all; close all; clc;


data = load('Fall.mat');


cells_idx = data.iscell(:,1) == 1;


F = data.F(cells_idx, :);     
Fneu = data.Fneu(cells_idx, :); 


F_corrected = F - 0.7 * Fneu;


F0 = prctile(F_corrected, 10, 2); 
F0_mat = repmat(F0, 1, size(F_corrected,2)); 


deltaF_over_F = (F_corrected - F0_mat) ./ F0_mat;


numFrames = size(F,2);
frameRate = 30; 
timeInSec = (0:numFrames-1) / frameRate;


figure('Position', [100, 100, 400, 400]);
imagesc(timeInSec, 1:size(deltaF_over_F,1), deltaF_over_F);
colormap('hot');
colorbar;


caxis([-0.05 0.5]);


xlabel('Time (s)');
ylabel('Neurons');
title('ΔF/F Neuronal Activity');
set(gca,'TickDir','out','FontSize',12);
axis tight;



