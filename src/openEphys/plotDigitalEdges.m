function fig = plotDigitalEdges(DIN, win_dur)
% PLOTDIGITALEDGES  Plot whole-session digital inputs as step traces built only from their
% transitions (DIN struct from importOpenEphysDigital or importIntanDigital), so no raw data is
% read. Inputs are stacked in one axes, first input on top. Only the visible window is drawn:
% where there are more transitions than pixels, each pixel-wide bin shows the range of states
% in it, so zooming in reveals every edge. A slider at the bottom scrolls a win_dur-s window
% through the session; pan/zoom also works. Opens on the full session; the button or slider
% switches to the window.
%
%   win_dur: visible window (s, default 30)

if ~exist('win_dur','var'); win_dur = 30; end
n_px = 2000;                                                     % bins across the visible window
n_chan = numel(DIN);
t_end = DIN(1).tend;

fig = figure('Name', 'Digital inputs');
ax = axes(fig); hold(ax, 'on')
labels = strings(n_chan, 1);
tr = cell(n_chan, 1); levels = cell(n_chan, 1); h = gobjects(n_chan, 1);
for k = 1:n_chan
    rise = DIN(k).data(:,1); rise = rise(~isnan(rise));
    fall = DIN(k).data(:,2); fall = fall(~isnan(fall));
    if isfield(DIN(1).info, 'initial_state'); start_high = DIN(1).info.initial_state(k);
    else; start_high = ~isempty(DIN(k).data) && isnan(DIN(k).data(1,1)); end
    tr{k} = sort([rise; fall]);
    levels{k} = double(xor(start_high, mod((0:numel(tr{k}))', 2) == 1)); % state before / after each transition
    h(k) = plot(ax, NaN, NaN, 'k');
    labels(n_chan - k + 1) = sprintf('%s: %d rises, %d falls', DIN(k).chanlabel, numel(rise), numel(fall));
end
set(ax, 'YTick', (0:n_chan-1) + 0.35, 'YTickLabel', labels, 'TickLabelInterpreter', 'none', ...
    'Box', 'off', 'TickDir', 'out', 'TickLength', [0.003 0.003])
ylim(ax, [-0.3 n_chan]); xlim(ax, [0 t_end]) % open on the full session
xlabel(ax, 'time from recording start (s)')

% Full session / window toggle (window position comes from the slider)
btn = uicontrol(fig, 'Style', 'togglebutton', 'Value', 1, 'String', sprintf('%g s window', win_dur), ...
    'Units', 'normalized', 'Position', [0.91 0.005 0.08 0.025], 'Callback', @toggleFull);
slider = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.2 0.005 0.7 0.02], ...
    'Min', 0, 'Max', t_end - win_dur, 'Value', 0, ...
    'SliderStep', [win_dur/2 win_dur] / (t_end - win_dur), 'Callback', @showWindow);
addlistener(slider, 'ContinuousValueChange', @showWindow); % update while dragging
ax.XAxis.LimitsChangedFcn = @(~,~) redraw();
redraw()
drawnow

    function redraw()
        lim = max(min(xlim(ax), t_end), 0);
        for kk = 1:n_chan
            row = n_chan - kk;                                   % first input on top
            t_k = tr{kk};
            i1 = find(t_k >= lim(1), 1); i2 = find(t_k <= lim(2), 1, 'last');
            state0 = levels{kk}(sum(t_k < lim(1)) + 1);          % state at the left edge
            if isempty(i1) || isempty(i2) || i2 < i1             % no transitions visible
                x = lim'; y = [state0; state0];
            elseif i2 - i1 + 1 <= n_px                           % exact steps
                x = [lim(1); repelem(t_k(i1:i2), 2, 1); lim(2)];
                y = repelem(levels{kk}(i1:i2+1), 2, 1);
            else                                                 % per-bin range of states
                edges = linspace(lim(1), lim(2), n_px + 1)';
                n_tr = histcounts(t_k(i1:i2), edges)';
                st = levels{kk}(discretize(edges(1:end-1), [-Inf; t_k; Inf]));  % state at each bin start
                y_lo = double(st == 0 | n_tr > 0); y_lo = 1 - y_lo;  % 0 if the bin reaches 0
                y_hi = double(st == 1 | n_tr > 0);                    % 1 if the bin reaches 1
                x = repelem(edges(1:end-1), 2);
                y = reshape([y_lo'; y_hi'], [], 1);
            end
            set(h(kk), 'XData', x, 'YData', y*0.7 + row)
        end
    end
    function showWindow(~, ~)
        xlim(ax, slider.Value + [0 win_dur]); btn.Value = 0; btn.String = 'Full session';
    end
    function toggleFull(~, ~)
        if btn.Value
            xlim(ax, [0 t_end]); btn.String = sprintf('%g s window', win_dur);
        else
            showWindow
        end
    end
end
