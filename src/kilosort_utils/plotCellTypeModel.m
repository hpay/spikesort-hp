function plotCellTypeModel(gm, res, plot_on, min_type_p, min_fit_p)
% PLOTCELLTYPEMODEL  Figures of the two-step whole-waveform model applied to units (res, from applyCellTypeModel
% 'waveform'). plot_on = 1: one column per cell type (E, I, nonsomatic, noise) if gm has cluster_type,
% otherwise one per cluster; 2: one column per cluster (titled with its cell type, if gm has them); 3: both.
% Click (without dragging) on a point in a 3D panel to highlight that unit in magenta: its waveform and ACG
% in its column, and its point in every 3D panel that shows it; its cluster id and session appear at the top
% left. Click and drag rotates.
% With an E / I split (gm.ei), the cell type view's 3D panels also mark E / I units whose E / I posterior is
% below min_type_p (orange x) and those whose fit to their E / I component has p below min_fit_p (black +):
% the cutoffs step 7 uses for good (both optional; omitted, or min_fit_p = 0, = no marks).
if ~exist('plot_on', 'var') || isempty(plot_on); plot_on = 1; end
if ~exist('min_type_p', 'var'); min_type_p = []; end
if ~exist('min_fit_p', 'var') || isempty(min_fit_p); min_fit_p = 0; end
types = {};
if isfield(res, 'types'); types = res.types; end   % cellTypesFromClusters (empty if the model has no cell types)
flag = struct('low', false(size(res.idx)), 'far', false(size(res.idx)));
if isfield(gm, 'ei') && ~isempty(types)
    is_ei = ismember(types(:), {'E', 'I'});
    flag.far = is_ei & res.ei_fit_p(:) < min_fit_p;
    if ~isempty(min_type_p); flag.low = is_ei & res.type_p(:) < min_type_p; end
end
flag.min_type_p = min_type_p; flag.min_fit_p = min_fit_p;
args = {res.idx, res.pcs, res.group, gm, res.basis, res.tWF, res.unit_labels, res.all_dir};
if plot_on ~= 1 || isempty(types); plotWaveformClusters(args{:}, {}, res.stats, res.z_label, res.ids, flag); end
if plot_on ~= 2 && ~isempty(types); plotWaveformClusters(args{:}, types, res.stats, res.z_label, res.ids, flag); end
end

function plotWaveformClusters(idx, pcs, group, gm, basis, tWF, unit_labels, all_dir, types, stats, z_label, ids, flag)
% One column per cell type (E, I, nonsomatic, noise) if types is given, otherwise per cluster:
% normalized peak-channel waveforms and autocorrelograms (good units thin black, MUA thin grey; waveform
% median in color; ACG means for good in color and MUA lighter).
% Bottom row, type view (4 panels; good units in the type's color, MUA lighter, open circles):
%   1. good E and I units: trough-to-peak width (ms), asymmetry, log10 rate (Hz)
%   2. good + MUA E and I units, same axes
%   3. E / I / nonsomatic units: wf PC1 (step 1 trough vs peak split), trough-group PC1 and PC3 (step 2)
%   4. same as 3 plus noise units: peak group black, trough group grey (noise column lines too)
% Bottom row, cluster view: each cluster's units (good in color, MUA lighter) among the other units of
% their wf PC1 group (grey) on that group's PCs 1 and 3 (x, y) and the z feature (default log rate), all
% z-scored; a column with units from both groups gets one small 3D panel per group.
% Bottom row: click and drag to rotate; linked views share a space.
is_good = strcmp(unit_labels(:), 'good');
type_view = exist('types', 'var') && ~isempty(types);
if type_view
    col_name = {'E', 'I', 'nonsomatic', 'noise'};
    col_mask = cellfun(@(c) strcmp(types(:), c), col_name, 'UniformOutput', false);
    colors = [0.85 0.1 0.1; 0.1 0.3 0.85; 0.1 0.6 0.2; 0.3 0.3 0.3];
    col_title = col_name;
    % Width, asymmetry (cellTypeFeatures) and log rate, in their own units
    ei_xyz = [basis.width_ms, basis.asym, basis.rate_log10];
    ei_lab = {'width (ms)', 'asymmetry', 'log_{10} rate (Hz)'};
    % Noise by group (noise cluster + outliers): peak group black, trough group grey
    is_noise = strcmp(types(:), 'noise'); is_outlier = is_noise & group(:) == 1;   % "outlier" = drawn grey
    w = basis.w_norm;
    m1 = gm.group(1); S1 = ((w - m1.mu) * m1.coeff - m1.s_mu) ./ m1.s_sd;  % every unit on the trough-group PCs
    model_xyz = [stats.wf_pc1, S1(:, [1 3])];
else
    n_clusters = max(idx); col_mask = arrayfun(@(k) idx == k, 1:n_clusters, 'UniformOutput', false);
    colors = lines(n_clusters);
    k_gauss = arrayfun(@(m) m.gm.NumComponents, gm.group);                  % background (if any) follows its group's Gaussians
    has_bg = arrayfun(@(m) isfield(m, 'p0') && m.p0 > 0, gm.group);
    has_out = arrayfun(@(m) isfield(m, 'nu') && ~isempty(m.nu) && isfield(m, 'outlier_p') && ~isempty(m.outlier_p) && m.outlier_p > 0, gm.group);
    extra = has_bg | has_out; col_title = compose('C%d', 1:n_clusters);
    extra_k = [k_gauss(1) + 1, k_gauss(1) + extra(1) + k_gauss(2) + 1]; extra_k = extra_k(extra);   % background / outlier clusters
    extra_name = {'background', 'outliers'}; extra_name = extra_name(1 + has_out(extra));
    if isfield(gm, 'cluster_type')   % final designation of each cluster, e.g. C1 (trough noise), C2 (excitatory)
        type_name = struct('E', 'excitatory', 'I', 'inhibitory', 'nonsomatic', 'nonsomatic', 'noise', 'noise');
        if isfield(gm, 'ei'); [type_name.E, type_name.I] = deal('somatic'); end   % E vs I then from width and asymmetry
        group_name = {'trough', 'peak'};
        for k = 1:n_clusters
            name = type_name.(gm.cluster_type{k});
            if strcmp(name, 'noise'); name = [group_name{1 + (k > k_gauss(1) + extra(1))} ' noise']; end
            if ismember(k, extra_k); name = [name ', ' extra_name{extra_k == k}]; end %#ok<AGROW>
            col_title{k} = sprintf('C%d (%s)', k, name);
        end
    else
        for i = 1:numel(extra_k); col_title{extra_k(i)} = sprintf('C%d (%s)', extra_k(i), extra_name{i}); end
    end
end
n_col = numel(col_mask);
t_ms = tWF * 1e3;
lag_ms = sqrt(basis.acg_edges(1:end-1) .* basis.acg_edges(2:end)) * 1e3;   % geometric bin centers
fig = figure('Name', 'Cell type model');
% Fixed axes positions (not tiledlayout, which re-fits every panel when a 3D panel is rotated)
col_w = 0.95 / n_col; row_y = [0.68 0.43 0.05]; row_h = [0.2 0.16 0.3];  % room for titles and x labels
pos = @(row, k) [0.045 + (k - 1)*col_w, row_y(row), 0.78*col_w, row_h(row)];
pc_ax = gobjects(0); pc_group = []; acg_ax = gobjects(0); wf_ax = gobjects(0);
for k = 1:n_col
    in_k = col_mask{k}; c = colors(k, :); c_light = c + 0.6*(1 - c);
    curves = {basis.w_norm, basis.acg};
    for row = 1:2
        ax = axes(fig, 'Position', pos(row, k)); hold(ax, 'on')
        x = {t_ms, lag_ms}; x = x{row};
        if type_view && strcmp(col_title{k}, 'noise')   % noise column by group: trough grey, peak black
            dark = in_k & ~is_outlier; light = in_k & is_outlier;
        else                                            % otherwise good black, MUA / noise grey
            dark = in_k & is_good; light = in_k & ~is_good;
        end
        if any(light); plot(ax, x, curves{row}(light, :)', 'Color', [0.6 0.6 0.6 0.4], 'LineWidth', 0.5); end
        if any(dark); plot(ax, x, curves{row}(dark, :)', 'Color', [0 0 0 0.4], 'LineWidth', 0.5); end
        if row == 1
            plot(ax, x, median(curves{row}(in_k, :), 1), 'Color', c, 'LineWidth', 2)
        else   % ACG: separate means for good (color) and MUA / noise (lighter), with a legend
            h_mean = gobjects(0); leg = {};
            if any(in_k & is_good); h_mean(end + 1) = plot(ax, x, mean(curves{row}(in_k & is_good, :), 1), 'Color', c, 'LineWidth', 2); leg{end + 1} = 'good'; end %#ok<AGROW>
            if any(in_k & ~is_good); h_mean(end + 1) = plot(ax, x, mean(curves{row}(in_k & ~is_good, :), 1), 'Color', c_light, 'LineWidth', 2); leg{end + 1} = 'mua'; end %#ok<AGROW>
            legend(ax, h_mean, leg, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8, 'AutoUpdate', 'off')
        end
        grid(ax, 'on'); xlim(ax, x([1 end]))
        if row == 1
            if type_view && strcmp(col_title{k}, 'noise')
                sub = {sprintf('%d peak (black), %d trough (grey)', nnz(in_k & ~is_outlier), nnz(in_k & is_outlier))};
            else   % short lines, so the axes toolbar at the top right doesn't cover them
                sub = {sprintf('%d good, %d mua', nnz(in_k & is_good), nnz(in_k & ~is_good))};
            end
            title(ax, [{sprintf('%s: n = %d', col_title{k}, nnz(in_k))}, sub], 'FontWeight', 'normal')
            if k == 1   % line colors for good and mua (stand-in lines, so the legend shows one of each)
                h_leg = [plot(ax, NaN, NaN, 'k'), plot(ax, NaN, NaN, 'Color', [0.6 0.6 0.6])];
                legend(ax, h_leg, {'good', 'mua'}, 'Location', 'southeast', 'Box', 'off', 'FontSize', 8, 'AutoUpdate', 'off')
            end
            ylim(ax, [-1.2 1.2]); xlabel(ax, 'time (ms)'); if k == 1; ylabel(ax, 'waveform / |largest extremum|'); end
            wf_ax(k) = ax; %#ok<AGROW>
        else
            set(ax, 'XScale', 'log'); xlabel(ax, 'lag (ms)'); acg_ax(k) = ax; %#ok<AGROW>
            if k == 1; ylabel(ax, sprintf('ACG: log_{10}(ratio to Poisson + %g)', basis.acg_offset)); end
        end
    end
    if type_view; continue; end     % type view: the bottom row is drawn below
    groups = unique(group(in_k))'; p3 = pos(3, k);
    for gi = 1:numel(groups)        % one 3D panel per group present (noise column: both)
        g = groups(gi); h = p3(4) / numel(groups);
        ax = axes(fig, 'Position', [p3(1), p3(2) + (numel(groups) - gi)*h, p3(3), h], 'PositionConstraint', 'innerposition');
        view(ax, 3); hold(ax, 'on') % 3D view before vis3d fixes the camera
        rest = group == g & ~in_k; mine = group == g & in_k;
        scatterUnits(ax, pcs, rest, 5, [0.7 0.7 0.7])
        scatterUnits(ax, pcs, mine & ~is_good, 8, c_light)
        scatterUnits(ax, pcs, mine & is_good, 8, c)
        ax.UserData = struct('xyz', pcs, 'shown', group == g);   % for highlighting a clicked unit
        group_name = {'trough cell predictors', 'peak cell predictors'};   % x, y = group PC1, PC3; z = z feature
        title(ax, group_name{g}, 'FontWeight', 'normal')
        xlabel(ax, 'PC1'); ylabel(ax, 'PC3'); zlabel(ax, [z_label{g} ' (z-scored)'])
        axis(ax, 'equal'); axis(ax, 'vis3d'); grid(ax, 'on')
        set(ax, 'XTickLabel', [], 'YTickLabel', [], 'ZTickLabel', [])           % no labels, so the panel can be bigger
        pc_ax(end + 1) = ax; pc_group(end + 1) = g; %#ok<AGROW>
    end
end
acg_shown = basis.acg(any([col_mask{:}], 2), :);
linkaxes(acg_ax, 'y'); ylim(acg_ax(1), prctile(acg_shown(:), [1 99]))   % shared ACG y range: 1st-99th percentile of all plotted values
if type_view
    is_ei = ismember(types(:), {'E', 'I'}); is_cell = ~strcmp(types(:), 'noise');
    model_lab = {'wf PC1 (split)', 'trough PC1', 'trough PC3'};
    panels = {ei_xyz, is_ei & is_good, false, 'good E/I', ei_lab; ...
              ei_xyz, is_ei, false, 'good+MUA E/I', ei_lab; ...
              model_xyz, is_cell, false, 'E/I/NS, model axes', model_lab; ...
              model_xyz, is_cell, true, '+ noise: peak (black), trough (grey)', model_lab};
    for k = 1:4
        [xyz, show, with_noise, ttl, lab] = panels{k, :};
        ax = axes(fig, 'Position', pos(3, k), 'PositionConstraint', 'innerposition'); view(ax, 3); hold(ax, 'on')
        xlabel(ax, lab{1}); ylabel(ax, lab{2}); zlabel(ax, lab{3})
        if with_noise   % noise by group: trough grey, peak black
            scatterUnits(ax, xyz, is_noise & is_outlier, 5, [0.6 0.6 0.6])
            scatterUnits(ax, xyz, is_noise & ~is_outlier, 5, [0 0 0])
            show = show | is_noise;
        end
        for t = 1:3   % E, I, nonsomatic: MUA lighter, then good in the type color
            in_t = show & col_mask{t}; c = colors(t, :);
            scatterUnits(ax, xyz, in_t & ~is_good, 8, c + 0.6*(1 - c))
            scatterUnits(ax, xyz, in_t & is_good, 8, c)
        end
        ax.UserData = struct('xyz', xyz, 'shown', show);         % for highlighting a clicked unit
        sel = show & flag.low;    % E / I posterior < min_type_p
        if any(sel); scatter3(ax, xyz(sel, 1), xyz(sel, 2), xyz(sel, 3), 40, [1 0.5 0], 'x', 'LineWidth', 1.5); end
        sel = show & flag.far;    % far from its E / I component
        if any(sel); scatter3(ax, xyz(sel, 1), xyz(sel, 2), xyz(sel, 3), 50, 'k', '+', 'LineWidth', 1.5); end
        title(ax, ttl, 'FontWeight', 'normal'); grid(ax, 'on')
        if k >= 3; axis(ax, 'equal'); end; axis(ax, 'vis3d')               % 1-2: different units on each axis
        set(ax, 'XTickLabel', [], 'YTickLabel', [], 'ZTickLabel', [])
        pc_ax(end + 1) = ax; pc_group(end + 1) = 1 + (k >= 3); %#ok<AGROW>
    end
end
for g = unique(pc_group)  % same viewing angle within the trough group and within the peak group
    setappdata(fig, sprintf('PCLink%d', g), linkprop(pc_ax(pc_group == g), 'View'))
end
set(findobj(fig, 'Type', 'axes'), 'FontSize', 10, 'TitleFontSizeMultiplier', 1, 'LabelFontSizeMultiplier', 1, 'Box', 'off', 'TickDir', 'out')
set([pc_ax.XLabel, pc_ax.YLabel, pc_ax.ZLabel], 'FontSize', 8)              % small axis labels on the 3D panels
rot = rotate3d(fig); rot.Enable = 'on';                                   % click and drag rotates the bottom row
setAllowAxesRotate(rot, setdiff(findobj(fig, 'Type', 'axes'), pc_ax), false)
% A click that doesn't rotate the view (no drag) selects the nearest unit
rot.ActionPreCallback = @(~, evd) setappdata(evd.Axes, 'view_at_press', evd.Axes.View);
rot.ActionPostCallback = @selectClickedUnit;
[parents, sessions] = fileparts(cellstr(all_dir));
[~, parents] = fileparts(parents);
label = annotation(fig, 'textbox', [0.005 0.955 0.4 0.04], 'String', 'click a point to highlight its unit', ...
    'EdgeColor', 'none', 'FontSize', 10, 'Color', 'm', 'Interpreter', 'none');
if type_view && any(flag.low | flag.far)
    key = {};
    if ~isempty(flag.min_type_p); key{end + 1} = sprintf('orange x: E/I posterior < %g (%d, %d good)', flag.min_type_p, nnz(flag.low), nnz(flag.low & is_good)); end
    if flag.min_fit_p > 0; key{end + 1} = sprintf('black +: fit to its E/I component p < %g (%d, %d good)', flag.min_fit_p, nnz(flag.far), nnz(flag.far & is_good)); end
    annotation(fig, 'textbox', [0.6 0.955 0.4 0.04], 'EdgeColor', 'none', 'FontSize', 10, 'HorizontalAlignment', 'right', ...
        'String', strjoin(key, '; '));
end
setappdata(fig, 'unit_info', struct('wf_ax', wf_ax, 'acg_ax', acg_ax, 'pc_ax', pc_ax, 'col_mask', {col_mask}, ...
    'col_title', {col_title}, 't_ms', t_ms, 'lag_ms', lag_ms, 'w_norm', basis.w_norm, 'acg', basis.acg, 'ids', ids, ...
    'session', {parents(basis.dir_idx)}, 'unit_labels', {unit_labels}, 'label', label))
model_name = 'Gaussian mixture'; if any(arrayfun(@(m) isfield(m, 'nu') && ~isempty(m.nu), gm.group)); model_name = 't-mixture'; end
n_good = nnz(is_good); n_mua = nnz(strcmp(unit_labels(:), 'mua'));
sgtitle(fig, sprintf('Two-step whole-waveform %s, %d clusters: %s %s (%d good + %d mua = %d)', model_name, max(idx), ...
    strjoin(unique(sessions), ', '), strjoin(unique(parents), ', '), n_good, n_mua, n_good + n_mua), 'Interpreter', 'none', 'FontSize', 10)
drawnow
end


function scatterUnits(ax, xyz, sel, sz, color)
% Open circles for the selected units
if ~any(sel); return; end
scatter3(ax, xyz(sel, 1), xyz(sel, 2), xyz(sel, 3), sz, color);
end

function selectClickedUnit(fig, evd)
% After a press in a 3D panel: if the view didn't change (a click, not a drag), highlight the unit whose
% point is nearest the line of sight through the click, if one is close (5% of the axes box)
ax = evd.Axes;
if ~isequal(ax.View, getappdata(ax, 'view_at_press')); return; end
d = ax.UserData; shown = find(d.shown);
lo = [ax.XLim(1), ax.YLim(1), ax.ZLim(1)]; scale = [diff(ax.XLim), diff(ax.YLim), diff(ax.ZLim)];
p = (d.xyz(shown, :) - lo) ./ scale;                   % points and click line in the unit axes box
a = (ax.CurrentPoint(1, :) - lo) ./ scale; b = (ax.CurrentPoint(2, :) - lo) ./ scale;
dist = vecnorm(cross(p - a, repmat(b - a, size(p, 1), 1), 2), 2, 2) / norm(b - a);
[d_min, i] = min(dist);
if d_min < 0.05; highlightUnit(fig, shown(i)); end
end

function highlightUnit(fig, u)
% Highlight unit u in magenta: waveform and ACG in its column, its point in every 3D panel that shows it
info = getappdata(fig, 'unit_info');
delete(findobj(fig, 'Tag', 'unit_highlight'))
k = find(cellfun(@(m) m(u), info.col_mask), 1);
plot(info.wf_ax(k), info.t_ms, info.w_norm(u, :), 'm', 'LineWidth', 2, 'Tag', 'unit_highlight')
plot(info.acg_ax(k), info.lag_ms, info.acg(u, :), 'm', 'LineWidth', 2, 'Tag', 'unit_highlight')
for a = info.pc_ax
    d = a.UserData;
    if d.shown(u); plot3(a, d.xyz(u, 1), d.xyz(u, 2), d.xyz(u, 3), 'mo', 'MarkerSize', 9, 'LineWidth', 2, 'Tag', 'unit_highlight'); end
end
info.label.String = sprintf('cluster %d, %s: %s, %s', info.ids(u), info.session{u}, info.col_title{k}, info.unit_labels{u});
end
