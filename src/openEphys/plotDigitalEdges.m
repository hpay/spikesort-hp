function fig = plotDigitalEdges(DIN, win_dur)
% PLOTDIGITALEDGES  Plot whole-session digital inputs as step traces built only from their
% transitions (DIN struct from importOpenEphysDigital or importIntanDigital), so no raw data is
% read. A slider at the bottom scrolls a win_dur-s window through the session; pan/zoom also
% works (x axes linked). Opens on the full session; the button or slider switches to the window.
%
%   win_dur: visible window (s, default 30)

if ~exist('win_dur','var'); win_dur = 30; end
n_chan = numel(DIN);
t_end = DIN(1).tend;

fig = figure('Name', 'Digital inputs', 'Position', [150 80 1400 110*n_chan]);
ax = gobjects(n_chan, 1);
for k = 1:n_chan
    rise = DIN(k).data(:,1); rise = rise(~isnan(rise));
    fall = DIN(k).data(:,2); fall = fall(~isnan(fall));
    if isfield(DIN(1).info, 'initial_state'); start_high = DIN(1).info.initial_state(k);
    else; start_high = ~isempty(DIN(k).data) && isnan(DIN(k).data(1,1)); end
    tr = sort([rise; fall]);
    levels = xor(start_high, mod((0:numel(tr))', 2) == 1); % state before / after each transition
    x = [0; repelem(tr, 2); t_end];
    y = repelem(double(levels), 2);
    ax(k) = subplot(n_chan, 1, k);
    plot(x, y, 'k')
    ylim([-0.2 1.2]); set(gca, 'YTick', [0 1], 'Box', 'off', 'TickDir', 'out', 'TickLength', [0.003 0.003])
    ylabel(DIN(k).chanlabel, 'Interpreter', 'none')
    title(sprintf('%s: %d rises, %d falls', DIN(k).chanlabel, numel(rise), numel(fall)), ...
        'FontSize', 8, 'Interpreter', 'none')
    if k < n_chan; set(gca, 'XTickLabel', []); end
end
xlabel('time from recording start (s)')
linkaxes(ax, 'x')
xlim(ax(1), [0 t_end]) % open on the full session

% Full session / window toggle (window position comes from the slider)
btn = uicontrol(fig, 'Style', 'togglebutton', 'Value', 1, 'String', sprintf('%g s window', win_dur), ...
    'Units', 'normalized', 'Position', [0.91 0.005 0.08 0.025], 'Callback', @toggleFull);
slider = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.13 0.005 0.775 0.02], ...
    'Min', 0, 'Max', t_end - win_dur, 'Value', 0, ...
    'SliderStep', [win_dur/2 win_dur] / (t_end - win_dur), 'Callback', @showWindow);
addlistener(slider, 'ContinuousValueChange', @showWindow); % update while dragging
drawnow

    function showWindow(~, ~)
        xlim(ax(1), slider.Value + [0 win_dur]); btn.Value = 0; btn.String = 'Full session';
    end
    function toggleFull(~, ~)
        if btn.Value
            xlim(ax(1), [0 t_end]); btn.String = sprintf('%g s window', win_dur);
        else
            showWindow
        end
    end
end
