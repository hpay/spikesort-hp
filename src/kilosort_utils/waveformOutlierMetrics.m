function [deriv_ss_out, thresh, wz] = waveformOutlierMetrics(mxWF, tWF, plot_on, n_sd, spike_win)
% WAVEFORMOUTLIERMETRICS  Noise metric on each unit's max-channel waveform (rows of mxWF),
% z-scored over time (as in the GMM_apply waveform plots).
%   deriv_ss_out: sum of squared sample-to-sample differences outside spike_win, so noisy,
%                 jagged waveforms score high but a sharp (narrow) spike itself doesn't count
%   thresh:       mean + n_sd*SD of log10(deriv_ss_out) (back in linear units)
% tWF: time (s) of each sample, trough at ~0
% plot_on: histogram (log x) plus waveforms of all units above thresh
% n_sd: default 2
% spike_win: [start stop] (s) around the trough excluded (default [-0.5 1] ms)

if ~exist('plot_on','var'); plot_on = false; end
if ~exist('n_sd','var'); n_sd = 2; end
if ~exist('spike_win','var'); spike_win = [-0.5 1]*1e-3; end

wz = zscore(mxWF, 0, 2);
t_mid = (tWF(1:end-1) + tWF(2:end))/2; % time of each difference
outside = t_mid < spike_win(1) | t_mid > spike_win(2);
dwz2 = diff(wz, 1, 2).^2;
deriv_ss_out = sum(dwz2(:, outside), 2);

log_vals = log10(deriv_ss_out(deriv_ss_out > 0));
thresh = 10^(mean(log_vals) + n_sd*std(log_vals));

if plot_on
    above = find(deriv_ss_out > thresh);
    [~, order] = sort(deriv_ss_out(above), 'descend');
    above = above(order);
    colors = orderedcolors('gem12'); colors = colors(mod(0:numel(above)-1, 12) + 1, :);

    figure('Position', [100 100 1400 550]);
    subplot(1, 2, 1)
    histogram(deriv_ss_out, logspace(min(log_vals), max(log_vals), 100), 'FaceColor', [.7 .7 .7]); hold on
    xline(thresh, '--k')
    for j = 1:numel(above) % mark each unit in its trace color
        plot(deriv_ss_out(above(j)), 0, 'v', 'MarkerFaceColor', colors(j,:), 'MarkerEdgeColor', 'none', 'MarkerSize', 7)
    end
    set(gca, 'XScale', 'log', 'XMinorTick', 'on')
    xlabel(sprintf('sum sq. derivatives outside [%g %g] ms (z)', spike_win*1e3)); ylabel('units')
    title(sprintf('n = %d; dashed = mean + %g SD of log10 = %.2f', nnz(~isnan(deriv_ss_out)), n_sd, thresh))

    subplot(1, 2, 2)
    h = gobjects(numel(above), 1);
    for j = 1:numel(above)
        h(j) = plot(tWF*1e3, wz(above(j),:), 'Color', colors(j,:)); hold on
    end
    plot(tWF*1e3, mean(wz, 1, 'omitnan'), 'k', 'LineWidth', 2)
    xline(spike_win*1e3, ':k')
    xlim(tWF([1 end])*1e3); ylim([-6 3])
    xlabel('time (ms)'); ylabel('z-scored amplitude')
    title(sprintf('%d units above threshold (black = mean of all units)', numel(above)))
    legend(h, compose('unit %d', above), 'FontSize', 7, 'Location', 'southeast')
    drawnow
end
end
