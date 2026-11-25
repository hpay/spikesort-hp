function fig = plotUnitBrowser(ks_dir, n_show, fig_title, groups)
% PLOTUNITBROWSER  Page through the good (or good + mua) units of a sort n_show at a time, to look for duplicates and
% bad splits. Units are ordered by shank pair (the two faces of one shank pair count as one x position),
% then by peak-channel depth, so units that could be the same neuron sit next to each other.
% Rows, per unit:
%   1. template on both faces of its shank pair within 100 um of its peak depth (thick = peak channel;
%      left = lower face number in channel_shanks.npy), over the next unit on the page in blue.
%      One amplitude scale for all units.
%   2. autocorrelogram, +/-50 ms, 1 ms bins
%   3. crosscorrelogram with the next unit in the order, and their template similarity
%      (similar_templates.npy, which Kilosort 4 uses to pick merge candidates)
% The slider's arrows and trough step n_show - 1 units, so the last unit of a page starts the next.
%
%   n_show (optional, default 8): units per page
%   fig_title (optional): e.g. session and sort name
%   groups (optional, default {'good'}): cluster_group.tsv labels to include, e.g. {'good', 'mua'}

if ~exist('n_show', 'var') || isempty(n_show); n_show = 8; end
if ~exist('fig_title', 'var') || isempty(fig_title); fig_title = ks_dir; end
if ~exist('groups', 'var'); groups = {'good'}; end
font_size = 8;
pair_pitch = 400;   % um: front faces are 400 um apart; each back face is within 400 um of its front

% Units, in shank pair then depth order
[good, wf, pos, peak_ch, amp_max, tmpl, group] = loadGoodTemplates(ks_dir, groups);
face = double(readNPY(fullfile(ks_dir, 'channel_shanks.npy')));
pair = floor(pos(:, 1) / pair_pitch);
[~, order] = sortrows([pair(peak_ch), pos(peak_ch, 2)]);
[good, wf, peak_ch, tmpl, group] = deal(good(order), wf(order, :, :), peak_ch(order), tmpl(order), group(order));
n_units = numel(good);

params = fileread(fullfile(ks_dir, 'params.py'));
fs = str2double(regexp(params, 'sample_rate\s*=\s*([\d.]+)', 'tokens', 'once'));
st = double(readNPY(fullfile(ks_dir, 'spike_times.npy'))) / fs;
clu = readNPY(fullfile(ks_dir, 'spike_clusters.npy'));
sim = readNPY(fullfile(ks_dir, 'similar_templates.npy'));
spk = arrayfun(@(c) st(clu == c), good, 'UniformOutput', false);
dur_s = max(st);

types = repmat({''}, n_units, 1);    % cell types from step 6 of runSingleKilosortOpenEphys.m, if run
if isfile(fullfile(ks_dir, 'gmm_result.mat'))
    r = load(fullfile(ks_dir, 'gmm_result.mat'));
    [found, loc] = ismember(good, r.id); types(found) = r.labels(loc(found));
end

wf_width = 20;                                  % um across for one waveform
face_gap = 50;                                  % um between the two faces in the drawing
amp_scale = 30 / prctile(amp_max, 90);          % um per template unit
tt = linspace(-wf_width/2, wf_width/2, size(wf, 2));
lags = -50:50;                                  % ms
next_color = [0 0.45 0.74];                     % next unit: template overlay and CCG

% Fixed axes (rows: templates, ACG, CCG)
fig = figure('Name', ['Unit browser: ' fig_title], 'Color', 'w');
row_y = [0.44 0.25 0.08]; row_h = [0.45 0.12 0.12];
col_w = 0.95 / n_show;
ax = gobjects(3, n_show);
for r = 1:3
    for c = 1:n_show
        ax(r, c) = axes(fig, 'Position', [0.045 + (c - 1)*col_w, row_y(r), 0.8*col_w, row_h(r)], ...
            'FontSize', font_size, 'TitleFontSizeMultiplier', 1, 'LabelFontSizeMultiplier', 1, ...
            'TickDir', 'out', 'Box', 'off', 'NextPlot', 'add');
    end
end
linkaxes(ax(1, :), 'y')                         % zooming one template panel zooms the row
set(ax(1, 2:end), 'YTickLabel', [])

n_first = max(n_units - n_show + 1, 1);         % last page starts here
sl = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.3 0.005 0.4 0.025], ...
    'Min', 1, 'Max', max(n_first, 2), 'Value', 1, ...
    'SliderStep', min((n_show - 1) * [1 1] / max(n_first - 1, 1), 1), 'Callback', @(s, ~) drawPage(round(s.Value)));
if n_first == 1; sl.Enable = 'off'; end
drawPage(1)

    function drawPage(first)
        sl.Value = first;
        arrayfun(@cla, ax); arrayfun(@(a) title(a, ''), ax);
        last = min(first + n_show - 1, n_units);
        y_drawn = [];
        for c = 1:n_show
            k = first + c - 1;
            set(ax(:, c), 'Visible', k <= n_units)
            if k > n_units; continue; end
            pk = peak_ch(k);

            % 1. templates on both faces of the shank pair
            in_pair = pair == pair(pk);
            faces = unique(face(in_pair));
            x_face = @(ch) pos(ch, 1) - min(pos(face == face(ch), 1)) + face_gap*(face(ch) ~= faces(1));
            near = find(in_pair & abs(pos(:, 2) - pos(pk, 2)) <= 100)';
            plot(ax(1, c), arrayfun(x_face, find(in_pair)), pos(in_pair, 2), '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 4)
            % next unit on this page (blue) underneath, on the same channels
            if k < last
                for ch = near
                    plot(ax(1, c), x_face(ch) + tt, pos(ch, 2) + squeeze(wf(k + 1, :, ch))*amp_scale, ...
                        'Color', next_color, 'LineWidth', 0.5)
                end
            end
            y_drawn = [y_drawn; pos(near, 2)]; %#ok<AGROW>
            for ch = near
                plot(ax(1, c), x_face(ch) + tt, pos(ch, 2) + squeeze(wf(k, :, ch))*amp_scale, 'k', ...
                    'LineWidth', 0.5 + (ch == pk))
            end
            set(ax(1, c), 'XLim', [-wf_width, face_gap + 16.5 + wf_width], 'XTick', [])
            title(ax(1, c), {sprintf('%d %s %s', good(k), types{k}, group{k}), ...
                sprintf('%.1f Hz, y %g', numel(spk{k})/dur_s, pos(pk, 2))}, 'FontWeight', 'normal', 'Interpreter', 'none')

            % 2. ACG, 3. CCG with the next unit
            K = ccg(spk{k}, spk{k}, 50, 1e-3); K(51) = 0;
            bar(ax(2, c), lags, K, 1, 'FaceColor', 'k', 'EdgeColor', 'none'); xlim(ax(2, c), lags([1 end]))
            if k < n_units
                K = ccg(spk{k}, spk{k + 1}, 50, 1e-3);
                bar(ax(3, c), lags, K, 1, 'FaceColor', next_color, 'EdgeColor', 'none'); xlim(ax(3, c), lags([1 end]))
                title(ax(3, c), sprintf('vs %d, sim %.2f', good(k + 1), sim(tmpl(k) + 1, tmpl(k + 1) + 1)), 'FontWeight', 'normal')
            end
        end
        set(ax(1, :), 'YLim', [min(y_drawn) - 15, max(y_drawn) + 15])  % channels drawn on this page
        ylabel(ax(1, 1), 'y (\mum)'); ylabel(ax(2, 1), 'ACG'); ylabel(ax(3, 1), 'CCG next'); xlabel(ax(3, 1), 'lag (ms)')
        sgtitle(fig, sprintf('%s: %s units %d-%d of %d (by shank pair, then depth)', fig_title, strjoin(groups, ' + '), ...
            first, last, n_units), 'FontSize', font_size, 'Interpreter', 'none')
        drawnow
    end
end
