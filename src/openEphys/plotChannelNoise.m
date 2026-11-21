function [noise_uV, bad_suggested] = plotChannelNoise(fbinary, nCh, fs, chan_map_file, bit_volts, n_chunks, snippet_dur)
% PLOTCHANNELNOISE  Per-channel high-pass noise from short snippets spread through
% a raw int16 binary (read only), plotted on the probe layout, to help pick bad channels.
% Also plots snippet_dur s of high-passed data from the middle of the session, channels stacked
% by shank and depth.
%
% noise_uV: SD per channel in the channel map, in map order (plain SD, so large repeated
%   artifacts raise it)
% bad_suggested: raw (1-indexed) channel numbers with noise > 2x or < 0.3x the median

if ~exist('n_chunks','var') || isempty(n_chunks); n_chunks = 30; end
if ~exist('snippet_dur','var'); snippet_dur = 4; end % s, from the middle of the recording
chunk_dur = 1; % s

m = load(chan_map_file);
f = dir(fbinary);
n_samp = f.bytes / (2*nCh);
[b, a] = butter(3, 300/(fs/2), 'high');

starts = round(linspace(fs*10, n_samp - fs*(chunk_dur+10), n_chunks));
mmf = memmapfile(fbinary, 'Format', {'int16', [nCh n_samp], 'x'});
sd_all = nan(numel(m.chanMap), n_chunks);
for ii = 1:n_chunks
    x = double(mmf.Data.x(m.chanMap, starts(ii) + (1:fs*chunk_dur)))' * bit_volts;
    sd_all(:,ii) = std(filtfilt(b, a, x));
end
noise_uV = median(sd_all, 2);
med_noise = median(noise_uV);
is_bad = noise_uV > 2*med_noise | noise_uV < 0.3*med_noise;
bad_suggested = sort(m.chanMap(is_bad))';

figure('Name', 'Channel noise', 'Position', [100 30 1400 1300]);
subplot(3,2,1)
scatter(m.xcoords, m.ycoords, 60, noise_uV, 'filled'); hold on
plot(m.xcoords(is_bad), m.ycoords(is_bad), 'rx', 'MarkerSize', 12, 'LineWidth', 2)
text(m.xcoords + 4, m.ycoords, string(m.chanMap), 'FontSize', 6)
c = colorbar; c.Label.String = 'noise SD (\muV, >300 Hz)';
xlabel('x (\mum)'); ylabel('y (\mum)'); title('Noise on probe layout (labels = raw CH#)')
subplot(3,2,2)
plot(m.chanMap, noise_uV, 'ko', 'MarkerFaceColor', 'k'); hold on
plot(m.chanMap(is_bad), noise_uV(is_bad), 'rx', 'MarkerSize', 12, 'LineWidth', 2)
yline([0.3 1 2]*med_noise, '--', {'0.3x median', 'median', '2x median'})
set(gca, 'YScale', 'log'); xlabel('raw channel #'); ylabel('noise SD (\muV)')
title(sprintf('Suggested bad: %s', mat2str(bad_suggested)))
grid on

% Raw snippet from the middle of the session, stacked by shank then depth (bad channels in red)
if snippet_dur > 0
    n_snip = round(fs*snippet_dur);
    snip_start = round(n_samp/2 - n_snip/2); % 0-indexed first sample
    x = filtfilt(b, a, double(mmf.Data.x(m.chanMap, snip_start + (1:n_snip)))' * bit_volts);
    [~, order] = sortrows([m.kcoords(:) m.ycoords(:)]);
    offset = 10*med_noise;
    t_snip = (snip_start + (0:n_snip-1))/fs;
    subplot(3,2,3:6); hold on
    for k = 1:numel(order)
        ch = order(k);
        plot(t_snip, x(:,ch) + k*offset, 'Color', [0.2 0.2 0.2] + 0.8*is_bad(ch)*[1 -0.2 -0.2], 'LineWidth', 0.25)
    end
    set(gca, 'YTick', (1:numel(order))*offset, 'YTickLabel', string(m.chanMap(order)), 'FontSize', 5)
    xlim(t_snip([1 end])); ylim([0 numel(order)+1]*offset)
    xlabel('time from recording start (s)'); ylabel('raw CH# (by shank, depth)')
    title(sprintf('%g s from mid-session, >300 Hz, spacing %.0f \\muV', snippet_dur, offset), 'FontSize', 10)
end
drawnow
