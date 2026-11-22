function [noise_uV, bad_suggested] = plotChannelNoise(fbinary, nCh, fs, chan_map_file, bit_volts, n_chunks, snippet_dur, session, zero_index)
% PLOTCHANNELNOISE  Per-channel high-pass noise from short snippets spread through
% a raw int16 binary (read only), to help pick bad channels. Plots snippet_dur s (default 30)
% from the middle of the session, channels stacked by shank and depth (see plotStackedTraces),
% with each channel's noise SD in a panel to the right, row for row.
%
% noise_uV: SD per channel in the channel map, in map order (plain SD, so large repeated
%   artifacts raise it)
% bad_suggested: raw (1-indexed) channel numbers with noise > 2x or < 0.5x the median
% session (optional): name shown at the start of the title
% zero_index (optional, default false): label channels 0-indexed (Intan / probe maker numbering)
%   instead of raw 1-indexed CH#, in all panels and the title (outputs stay 1-indexed)

if ~exist('n_chunks','var') || isempty(n_chunks); n_chunks = 30; end
if ~exist('snippet_dur','var') || isempty(snippet_dur); snippet_dur = 30; end % s, from the middle of the recording
if ~exist('session','var'); session = ''; end
if ~exist('zero_index','var'); zero_index = false; end
ch_offset = double(zero_index); % subtracted from raw CH# in labels
ch_name = 'raw CH#'; if zero_index; ch_name = 'ch (0-indexed)'; end
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
is_bad = noise_uV > 2*med_noise | noise_uV < 0.5*med_noise;
bad_suggested = sort(m.chanMap(is_bad))';

% Channel order: by shank, then depth. Double-sided probe: front and back adjacent, in the map's
% virtual shank order if it already alternates front/back, otherwise A_front, A_back, B_front, ...
[~, first, group] = unique(m.kcoords(:));
group_names = "shank " + string(m.kcoords(first));
if isfield(m, 'shank_letter')
    names = string(m.shank_letter) + "_" + replace(string(m.side), ["f" "b"], ["front" "back"]);
    if ~all(abs(diff(m.side(first) == 'b')) == 1)
        [~, first, group] = unique(2*double(m.shank_letter) + (m.side == 'b'));
    end
    group_names = names(first);
end
[~, order] = sortrows([group m.ycoords(:)]);
rows = m.chanMap(order);
n_ch = numel(order);

font_size = 10;  % all text except the per-channel labels
label_size = 6;  % per-channel labels (trace rows)
tick_pt = 3;     % tick length (points) on every axes
fig = figure('Name', 'Channel noise');
sgtitle(fig, sprintf('%s: Suggested bad %s: %s', session, ch_name, mat2str(bad_suggested - ch_offset)), 'Interpreter', 'none', 'FontSize', font_size + 2) % plotStackedTraces appends the trace display mode

% Traces, opening on snippet_dur s from the middle of the session. Interactive: slider through
% the whole recording, high-pass / wide band, median subtraction and hide-bad-channel toggles
% (see plotStackedTraces).
if snippet_dur > 0
    read_fun = @(i0, n) readRows(fbinary, nCh, rows, i0, n) * bit_volts;
    ax = axes(fig, 'Position', [0.05 0.06 0.69 0.83]);
    plotStackedTraces(ax, read_fun, fs, n_samp, round(n_samp/2)/fs - snippet_dur/2, snippet_dur, ...
        string(rows - ch_offset), repmat(0.2, n_ch, 3), group(order), is_bad(order))
    ylabel(ax, [ch_name ' (by shank, depth)'])
end

% Noise SD per channel, same rows as the traces
ax_n = axes(fig, 'Position', [0.755 0.06 0.065 0.83]); hold(ax_n, 'on')
plot(ax_n, noise_uV(order), 1:n_ch, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 3)
bad_row = find(is_bad(order));
plot(ax_n, noise_uV(order(bad_row)), bad_row, 'rx', 'MarkerSize', 10, 'LineWidth', 2)
xline(ax_n, [0.5 1 2]*med_noise, '--')
yline(ax_n, find(diff(group(order)) ~= 0) + 0.5, 'Color', 'c', 'LineWidth', 1.5)
set(ax_n, 'XScale', 'log', 'YTick', [], 'Box', 'off', 'TickDir', 'out')
xlim(ax_n, [min(0.45*med_noise, 0.9*min(noise_uV)), max(2.2*med_noise, 1.1*max(noise_uV))]) % just below the 0.5x cutoff, unless a channel is lower
ylim(ax_n, [0 n_ch+1]); grid(ax_n, 'on')
xlabel(ax_n, {'noise SD (\muV, >300 Hz)', 'lines: 0.5, 1, 2x median'})
y_mid = accumarray(group(order), (1:n_ch)', [], @mean); % middle row of each group
text(ax_n, ones(size(y_mid)), y_mid/(n_ch + 1), "  " + group_names, 'Units', 'normalized', ...
    'Interpreter', 'none', 'FontSize', font_size, 'FontWeight', 'bold')
if snippet_dur > 0; linkaxes([ax ax_n], 'y'); end

% Contact locations from the channel map, rotated (map x up, depth y across) so shanks stack
% like the traces; labels = channel numbers, bad channels in red
ax_p = axes(fig, 'Position', [0.875 0.06 0.115 0.83]); hold(ax_p, 'on')
plot(ax_p, m.ycoords, m.xcoords, 's', 'MarkerSize', 13, 'Color', [0.8 0.8 0.8], 'MarkerFaceColor', [0.9 0.9 0.9])
t = text(ax_p, m.ycoords, m.xcoords, string(m.chanMap - ch_offset), 'FontSize', font_size - 1, 'HorizontalAlignment', 'center');
set(t(is_bad), 'Color', 'r')
xlim(ax_p, [min(m.ycoords) - 10, max(m.ycoords) + 10]); ylim(ax_p, [min(m.xcoords) - 30, max(m.xcoords) + 40])
set(ax_p, 'YTick', [], 'Box', 'off', 'TickDir', 'out', 'YColor', 'none')
xlabel(ax_p, 'map y (\mum)'); title(ax_p, ['contacts (' ch_name ')'], 'FontWeight', 'normal')

% Same font size and tick length (in points) everywhere
drawnow % let the (docked) figure reach its final size first
for a = findobj(fig, 'Type', 'axes')'
    a.FontSize = font_size;
    units = a.Units; a.Units = 'points'; pos = a.Position; a.Units = units;
    a.TickLength = [1 1] * tick_pt / max(pos(3:4));
end
if snippet_dur > 0; ax.YAxis.FontSize = label_size; ax.YLabel.FontSize = font_size; end
drawnow
end

function x = readRows(fbinary, nCh, rows, i0, n)
% n samples of the given raw channels (rows, 1-indexed) starting at 0-indexed sample i0,
% as samples x channels, int16 converted to single
fid = fopen(fbinary, 'r');
fseek(fid, 2*nCh*i0, 'bof');
x = fread(fid, [nCh n], '*int16');
fclose(fid);
x = single(x(rows, :))';
end
