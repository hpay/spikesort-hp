function fig = plotDriftmap(ks_dir, fs, t_bin, y_bin)
% PLOTDRIFTMAP  Spike density vs. time and depth over the whole session, to check for drift.
% Works the same way for Kilosort 2 and 4 output. Each spike's PC-feature energy (sum of
% squares over PCs, so sign doesn't matter) is computed on each of its template's channels;
% its depth is the peak of a Gaussian fit through the amplitude (sqrt energy) on its largest
% channel and the channels one row above and below on the same shank. Shanks are stacked bottom to top;
% within a shank, depth increases up. The image is log(1 + spike count) in t_bin (s, default 10)
% x y_bin (um, default 2) bins.

if ~exist('t_bin','var'); t_bin = 10; end
if ~exist('y_bin','var'); y_bin = 2; end

st = double(readNPY(fullfile(ks_dir, 'spike_times.npy'))) / fs;
pos = readNPY(fullfile(ks_dir, 'channel_positions.npy'));       % channels x [x y], channel_map order
spike_templates = readNPY(fullfile(ks_dir, 'spike_templates.npy'));
pc_ind = readNPY(fullfile(ks_dir, 'pc_feature_ind.npy')) + 1;   % templates x channels
pc = readNPY(fullfile(ks_dir, 'pc_features.npy'));              % spikes x PCs x channels
energy = squeeze(sum(single(pc).^2, 2)); clear pc               % spikes x channels
chans = pc_ind(double(spike_templates) + 1, :);

% Shank of each channel: x positions more than 50 um apart are different shanks
[xs, ~, ix] = unique(pos(:,1));
shank_of_x = cumsum([1; diff(xs) > 50]);
shank = shank_of_x(ix);
row_pitch = median(diff(unique(pos(:,2))));                      % vertical channel spacing
y_ch = single(pos(:,2));

% Sub-row peak position: fit a Gaussian through the amplitude on the largest channel and the
% channels one row above and below it on the same shank (3-point interpolation). Unlike a
% center of mass, this doesn't pile spikes up at contact depths.
amp = sqrt(energy); clear energy
[a0, k_max] = max(amp, [], 2);
max_chan = chans(sub2ind(size(chans), (1:numel(k_max))', k_max));
same_shank = shank(chans) == shank(max_chan);
dy = y_ch(chans) - y_ch(max_chan);
a_up = sum(amp .* (same_shank & abs(dy - row_pitch) < 1), 2);   % 0 if that channel isn't among the template's
a_dn = sum(amp .* (same_shank & abs(dy + row_pitch) < 1), 2);
clear amp same_shank dy
offset = zeros(size(a0), 'single');
both = a_up > 0 & a_dn > 0;
l0 = log(a0(both)); lu = log(a_up(both)); ld = log(a_dn(both));
offset(both) = row_pitch * (lu - ld) ./ (2 * (2*l0 - lu - ld));
offset = max(min(offset, row_pitch/2), -row_pitch/2);
depth = y_ch(max_chan) + offset;

y_span = max(pos(:,2)) - min(pos(:,2)) + y_bin;                  % shank height, shanks stacked with no gap
y_plot = (shank(max_chan) - 1) * y_span + double(depth);

y_edges = min(pos(:,2)) : y_bin : max(shank)*y_span + min(pos(:,2));
counts = histcounts2(st, y_plot, 0:t_bin:max(st)+t_bin, y_edges);
fig = figure('Name', 'Drift map');
imagesc([0 max(st)], y_edges([1 end-1]) + y_bin/2, log1p(counts)')
axis xy; colormap(flipud(gray)); c = colorbar; c.Label.String = sprintf('log(1 + spikes per %g s x %g \\mum)', t_bin, y_bin);
n_shank = max(shank);
set(gca, 'YTick', ((1:n_shank) - 0.5) * y_span + min(pos(:,2)), 'YTickLabel', compose('shank %d', 1:n_shank), ...
    'TickDir', 'out', 'Box', 'off')
hold on; yline((1:n_shank-1) * y_span + min(pos(:,2)), 'Color', [0.8 0.3 0.3])
xlabel('time (s)'); ylabel('depth within shank (bottom to top)')
[~, sort_name] = fileparts(ks_dir); [~, session] = fileparts(fileparts(ks_dir));
title(sprintf('%s %s: spike depth (3-point Gaussian peak) over time', session, sort_name), 'Interpreter', 'none')
drawnow
end
