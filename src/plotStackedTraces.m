function plotStackedTraces(ax, read_fun, fs, n_samp, t_start, win_dur, labels, colors, group, bad)
% PLOTSTACKEDTRACES  Interactive stacked multichannel traces from a long recording, read on
% demand. Only the visible time window is read (plus a margin so filtering has no edge
% effects) and drawn, reduced to a min/max envelope of ~2000 points per channel; zooming or
% panning reads and redraws as needed. A slider at the top right moves through the whole
% recording; its arrows jump by the current x-range and clicking its track by 5x that.
% Buttons at the top left toggle
%   - high-pass (300 Hz, 3rd-order Butterworth, zero phase) vs. wide band (default wide band)
%   - subtracting the median across channels at each sample (common reference; default on)
% The trace spacing is 10x the median channel SD, set when the display mode first changes.
%
%   ax: axes to draw in
%   read_fun(i0, n): n samples starting at 0-indexed sample i0, as samples x channels (uV),
%       columns in display order (bottom to top)
%   fs: sample rate; n_samp: samples in the recording
%   t_start, win_dur: initial window (s); labels: y tick label per channel; colors: channels x 3 RGB
%   group (optional): group of each channel (e.g. shank); cyan lines are drawn between groups
%   bad (optional): logical per channel; adds a button that hides these channels' traces

n_px = 2000;                                     % points per channel to draw
max_dur = 60;                                    % longest window read at once (s)
pad = round(0.1*fs);                             % read margin for the filter (samples)
t_end = (n_samp - 1)/fs;
[b, a] = butter(3, 300/(fs/2), 'high');
state = struct('high', false, 'car', true);
spacing = containers.Map();                      % trace spacing for each display mode
loaded = [];                                     % sample range of the data in memory
y_all = [];                                      % displayed signal for that range
n_ch = numel(labels);

axes(ax); cla(ax); hold(ax, 'on')
h = gobjects(n_ch, 1);
for k = 1:n_ch
    h(k) = plot(ax, NaN, NaN, 'Color', colors(k,:), 'LineWidth', 0.5);
end
set(ax, 'YTick', 1:n_ch, 'YTickLabel', labels, 'TickDir', 'out', 'Box', 'off')
ylim(ax, [0 n_ch+1])
if exist('group', 'var')
    yline(ax, find(diff(group(:)) ~= 0) + 0.5, 'Color', 'c', 'LineWidth', 1.5)
end
xlabel(ax, 'time from recording start (s)')

fig = ancestor(ax, 'figure'); p = ax.Position; y_btn = p(2) + p(4) + 0.005;
sg = findall(fig, 'Type', 'subplottext'); % if the figure has an sgtitle, the display mode is appended to it
if ~isempty(sg); sg_base = string(sg.String); end
uicontrol(fig, 'Style', 'togglebutton', 'String', 'High-pass 300 Hz', 'Value', 0, 'Units', 'normalized', ...
    'Position', [p(1) y_btn 0.1 0.03], 'Callback', @(s,~) setState('high', s.Value));
uicontrol(fig, 'Style', 'togglebutton', 'String', 'Subtract median', 'Value', 1, 'Units', 'normalized', ...
    'Position', [p(1)+0.105 y_btn 0.1 0.03], 'Callback', @(s,~) setState('car', s.Value));
if exist('bad', 'var') && any(bad)
    uicontrol(fig, 'Style', 'togglebutton', 'String', 'Hide bad chans', 'Value', 0, 'Units', 'normalized', ...
        'Position', [p(1)+0.21 y_btn 0.1 0.03], 'Callback', @(s,~) set(h(bad), 'Visible', ~s.Value));
end
% Slider over the whole recording: arrows jump by the current x-range, the track by 5x
slider = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', 'Position', [p(1)+p(3)-0.35 y_btn+0.005 0.25 0.02], ... % right corner left free for the axes toolbar
 ...
    'Min', 0, 'Max', t_end, 'Value', t_start, 'Callback', @(s,~) moveTo(s.Value));

xlim(ax, t_start + [0 win_dur]);
ax.XAxis.LimitsChangedFcn = @(~,~) redraw();
redraw()

    function setState(name, value)
        state.(name) = logical(value);
        loaded = [];                             % reprocess for the new mode
        redraw()
    end

    function moveTo(t1)
        w = diff(xlim(ax));
        t1 = min(max(t1, 0), max(t_end - w, 0));
        xlim(ax, t1 + [0 w]);                    % triggers redraw
    end

    function load(i1, i2)
        % read samples i1:i2 (0-indexed) with a margin, filter, subtract the median, trim
        r0 = max(i1 - pad, 0); r1 = min(i2 + pad, n_samp - 1);
        x = single(read_fun(r0, r1 - r0 + 1));
        x = x - mean(x);
        if state.high; x = filtfilt(b, a, x); end
        if state.car; x = x - median(x, 2); end
        y_all = x((i1 - r0) + 1 : (i2 - r0) + 1, :);
        loaded = [i1 i2];
    end

    function redraw()
        lim = xlim(ax);
        if diff(lim) > max_dur                   % don't read more than max_dur at once
            lim(2) = lim(1) + max_dur;
            xlim(ax, lim); return                % redraws through LimitsChangedFcn
        end
        i1 = max(0, floor(lim(1)*fs)); i2 = min(n_samp - 1, ceil(lim(2)*fs));
        if isempty(loaded) || i1 < loaded(1) || i2 > loaded(2)
            load(i1, i2);
        end
        key = sprintf('%d%d', state.high, state.car);
        if ~isKey(spacing, key); spacing(key) = 10 * median(std(y_all)); end
        y = y_all((i1:i2) - loaded(1) + 1, :);
        tt = (i1:i2)'/fs;
        if size(y, 1) > 2*n_px                   % min/max envelope per block
            blk = floor(size(y, 1) / n_px);
            nb = floor(size(y, 1) / blk);
            yb = reshape(y(1:nb*blk, :), blk, nb, n_ch);
            y = reshape([min(yb, [], 1); max(yb, [], 1)], 2*nb, n_ch);
            tb = tt(1:blk:nb*blk);
            tt = reshape([tb'; tb'], [], 1);
        end
        for kk = 1:n_ch
            set(h(kk), 'XData', tt, 'YData', y(:, kk) / spacing(key) + kk)
        end
        step = min(diff(lim) / t_end, 1);       % arrows: one x-range; track: 5x
        set(slider, 'Value', min(max(lim(1), 0), t_end), 'SliderStep', [step min(5*step, 1)]);
        mode_str = 'wide band'; if state.high; mode_str = '>300 Hz'; end
        if state.car; mode_str = [mode_str ', median subtracted']; end
        mode_str = sprintf('%s, spacing %.0f %sV', mode_str, spacing(key), char(181));
        if isempty(sg); title(ax, mode_str); else; sg.String = sg_base + "   (traces: " + mode_str + ")"; end
    end
end
