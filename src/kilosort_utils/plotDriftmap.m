function fig = plotDriftmap(ks_dir, fs, t_bin, y_bin)
% PLOTDRIFTMAP  Spike density vs. time and probe position over the whole session, to check
% for drift. Works for Kilosort 2 and 4 output. Each spike's depth is the center of mass of
% its first PC feature over its template's channels on the template's shank (as in the
% cortex-lab ksDriftmap). Shanks are stacked bottom to top; within a shank, depth increases up.
% The image is log(1 + spike count) in t_bin (s, default 10) x y_bin (um, default 2) bins.

if ~exist('t_bin','var'); t_bin = 10; end
if ~exist('y_bin','var'); y_bin = 2; end

st = double(readNPY(fullfile(ks_dir, 'spike_times.npy'))) / fs;
spike_templates = readNPY(fullfile(ks_dir, 'spike_templates.npy'));
pos = readNPY(fullfile(ks_dir, 'channel_positions.npy'));       % channels x [x y], channel_map order
pc_ind = readNPY(fullfile(ks_dir, 'pc_feature_ind.npy')) + 1;   % templates x channels, first = main channel
pc = readNPY(fullfile(ks_dir, 'pc_features.npy'));              % spikes x PCs x channels
pc1 = max(single(squeeze(pc(:,1,:))), 0); clear pc              % first PC only, negatives ignored

% Shank of each channel: x positions more than 50 um apart are different shanks
[xs, ~, ix] = unique(pos(:,1));
shank_of_x = cumsum([1; diff(xs) > 50]);
shank = shank_of_x(ix);

chans = pc_ind(double(spike_templates) + 1, :);                 % spikes x channels
w = pc1.^2 .* (shank(chans) == shank(chans(:,1)));              % only channels on the main channel's shank
y_ch = pos(:,2);
depth = sum(y_ch(chans) .* w, 2) ./ sum(w, 2);
y_span = max(pos(:,2)) - min(pos(:,2)) + 50;                     % shank height + gap
y_plot = (shank(chans(:,1)) - 1) * y_span + depth;

counts = histcounts2(st, y_plot, 0:t_bin:max(st)+t_bin, min(y_plot):y_bin:max(y_plot)+y_bin);
fig = figure('Name', 'Drift map', 'Position', [100 100 1400 800]);
imagesc([0 max(st)], [min(y_plot) max(y_plot)], log1p(counts)')
axis xy; colormap(flipud(gray)); c = colorbar; c.Label.String = sprintf('log(1 + spikes per %g s x %g \\mum)', t_bin, y_bin);
n_shank = max(shank);
set(gca, 'YTick', ((1:n_shank) - 0.5) * y_span + min(pos(:,2)) - 25, 'YTickLabel', compose('shank %d', 1:n_shank), ...
    'TickDir', 'out', 'Box', 'off')
hold on; yline((1:n_shank-1) * y_span + min(pos(:,2)) - 25, 'Color', [0.8 0.3 0.3])
xlabel('time (s)'); ylabel('depth within shank (bottom to top)')
[~, sort_name] = fileparts(ks_dir); [~, session] = fileparts(fileparts(ks_dir));
title(sprintf('%s %s: spike positions over time', session, sort_name), 'Interpreter', 'none')
drawnow
end
