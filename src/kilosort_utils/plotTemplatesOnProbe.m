function fig = plotTemplatesOnProbe(ks_dir, fig_title)
% PLOTTEMPLATESONPROBE  Template waveforms of the final good units (cluster_group.tsv), drawn at the
% channel positions the sort used (channel_positions.npy), one panel per shank. Shanks are channel
% groups more than 50 um apart in x, so faces placed at the same x (aligned map) share a panel.
% Within a panel each unit has its own lane (ordered by peak depth, cluster id on top), and its
% waveforms on the channels within 50 um of its peak channel are drawn at their true position
% within the shank (x relative to the shank's leftmost column, y). Templates are unwhitened
% (whitening_mat_inv.npy); one amplitude scale for all units. Works for Kilosort 2 and 4 output.
%
%   fig_title (optional): figure title, e.g. session and sort name

[good, wf, pos, peak_ch, amp_max] = loadGoodTemplates(ks_dir);
[n_samp, n_ch] = deal(size(wf, 2), size(wf, 3));

% Shanks: channel groups separated by more than 50 um in x
[x_sorted, ix] = sort(pos(:, 1));
shank = zeros(n_ch, 1); shank(ix) = cumsum([1; diff(x_sorted) > 50]);
n_shank = max(shank);

wf_width = 14;                                  % um across for one waveform
amp_scale = 20 / prctile(amp_max, 90);          % um per template unit (90th percentile peak spans 20 um)
tt = linspace(-wf_width/2, wf_width/2, n_samp);
colors = lines(7);

fig = figure('Name', 'Templates on probe');
tl = tiledlayout(fig, 1, n_shank, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(n_shank, 1);
for s = 1:n_shank
    ax(s) = nexttile(tl); hold(ax(s), 'on')
    on_shank = shank == s;
    x0 = min(pos(on_shank, 1));
    lane_w = max(pos(on_shank, 1)) - x0 + wf_width + 10;          % one unit per lane
    [~, order] = sort(pos(peak_ch, 2));
    units = order(shank(peak_ch(order)) == s)';                   % this shank's units, by peak depth
    for li = 1:numel(units)
        k = units(li);
        near = find(on_shank & vecnorm(pos - pos(peak_ch(k), :), 2, 2) <= 50)';
        c = colors(mod(li - 1, 7) + 1, :);
        xl = (li - 1)*lane_w + pos(on_shank, 1) - x0;
        plot(ax(s), xl, pos(on_shank, 2), '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 4)
        xx = (li - 1)*lane_w + pos(near, 1)' - x0 + tt';
        yy = pos(near, 2)' + squeeze(wf(k, :, near)) * amp_scale;
        plot(ax(s), xx, yy, 'Color', c, 'LineWidth', 0.75)
        text(ax(s), (li - 1)*lane_w + (lane_w - 10)/2, max(pos(:, 2)) + 20, sprintf('%d', good(k)), ...
            'Color', c, 'FontSize', 7, 'HorizontalAlignment', 'center')
    end
    xlim(ax(s), [-wf_width, max(numel(units), 1)*lane_w])
    title(ax(s), sprintf('x %g um (n=%d)', x0, numel(units)), 'FontWeight', 'normal', 'FontSize', 9)
    set(ax(s), 'TickDir', 'out', 'Box', 'off', 'FontSize', 9, 'XTick', [])
    if s > 1; ax(s).YTickLabel = []; end
end
linkaxes(ax, 'y'); ylim(ax(1), [min(pos(:, 2)) - 25, max(pos(:, 2)) + 25])
ylabel(ax(1), 'y (\mum)'); xlabel(tl, 'units on each shank (lanes ordered by peak depth); panel title = shank x')
if ~exist('fig_title', 'var'); fig_title = ks_dir; end
title(tl, sprintf('%s: templates of %d good units', fig_title, numel(good)), 'Interpreter', 'none')
drawnow
end
