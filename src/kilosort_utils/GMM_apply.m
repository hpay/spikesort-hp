function [idx, labels, stats, goodIDs, P] = GMM_apply(gm, all_dir, option_only_good, plot_on, stat_set, basis)
% stat_set: 'trough' (default) or 'waveform' (see GMM_stats); basis: from GMM_make ('waveform')
% 'waveform': gm is the two-step struct from GMM_make. idx = cluster number, numbered group 1 (negative
%   trough) first, then group 2 (0 if features missing); labels 'C1', 'C2', ..., or the cell type
%   (E / I / nonsomatic / noise; missing features = noise) if the model has a cluster_type map, and then
%   the figure shows one column per cell type

if ~exist('option_only_good','var'); option_only_good = []; end
if ~exist('stat_set', 'var') || isempty(stat_set); stat_set = 'trough'; end
if ~exist('basis', 'var'); basis = []; end
id = 'stats:gmdistribution:cluster:MissingData'; warning('off',id)

if strcmp(stat_set, 'waveform')
    [stats, goodIDs, ~, tWF, ~, unit_labels, basis] = GMM_stats(all_dir, option_only_good, stat_set, basis);
    group = 2 - (gm.split.cluster(stats.wf_pc1) == gm.split_low);
    X_other = stats{:, gm.other};
    has_bg = arrayfun(@(m) isfield(m, 'p0') && m.p0 > 0, gm.group);
    % t models: outlier_p > 0 moves units whose fit to their own component has p < outlier_p (squared
    % distance / d ~ F(d, nu)) to an outlier cluster, which takes the background's place
    has_out = arrayfun(@(m) isfield(m, 'nu') && ~isempty(m.nu) && isfield(m, 'outlier_p') && ~isempty(m.outlier_p) && m.outlier_p > 0, gm.group);
    k_group = arrayfun(@(m) m.gm.NumComponents, gm.group) + (has_bg | has_out);   % background / outliers: the group's last cluster
    idx = zeros(height(stats), 1); P = zeros(height(stats), sum(k_group)); pcs = NaN(height(stats), 3);
    for g = 1:2
        m = gm.group(g); in_g = group == g;
        pc_g = ((basis.w_norm(in_g, :) - m.mu) * m.coeff - m.s_mu) ./ m.s_sd;   % group PCs 1-3, z-scored
        F = [pc_g(:, 1:m.n_pc), X_other(in_g, :)];
        nu = []; if isfield(m, 'nu'); nu = m.nu; end                   % nu: t-distribution components
        if has_bg(g) || ~isempty(nu)   % posterior over the components (and the uniform background)
            L = (1 - m.p0) * m.gm.ComponentProportion .* componentDensity(F, m.gm, nu);
            if has_bg(g); L = [L, repmat(m.p0 * m.bg_density, size(F, 1), 1)]; end %#ok<AGROW>
            P_g = L ./ sum(L, 2); [~, idx_g] = max(P_g, [], 2); idx_g(any(isnan(F), 2)) = NaN;
            if has_out(g)
                K = m.gm.NumComponents; d = size(F, 2); p_fit = NaN(size(idx_g));
                for k = 1:K
                    in_k = idx_g == k; Xc = F(in_k, :) - m.gm.mu(k, :);
                    p_fit(in_k) = 1 - fcdf(sum((Xc / m.gm.Sigma(:, :, min(k, end))) .* Xc, 2) / d, d, nu);
                end
                idx_g(p_fit < m.outlier_p) = K + 1; P_g(:, K + 1) = 0;
            end
        else
            [idx_g, ~, P_g] = m.gm.cluster(F);
        end
        if isfield(m, 'plot_z') && ~isempty(m.plot_z); z_name = m.plot_z; else; z_name = 'rate_log10'; end
        pcs(in_g, :) = [pc_g(:, [1 3]), stats.(z_name)(in_g)];          % plotted: PC1, PC3, z feature (default log rate)
        offset = sum(k_group(1:g-1));
        idx(in_g) = idx_g + offset; P(in_g, offset + (1:k_group(g))) = P_g;
    end
    warning('on', id)
    idx(isnan(idx)) = 0;
    labels = compose('C%d', idx); labels(idx == 0) = {'unknown'};
    for k = 1:sum(k_group)
        fprintf('C%d: %d units (%d good, %d mua, %d noise)\n', k, nnz(idx == k), nnz(idx == k & strcmp(unit_labels, 'good')), ...
            nnz(idx == k & strcmp(unit_labels, 'mua')), nnz(idx == k & strcmp(unit_labels, 'noise')))
    end
    types = {};
    if isfield(gm, 'cluster_type')   % model maps clusters to cell types: labels are E / I / nonsomatic / noise
        types = repmat({'noise'}, size(idx)); types(idx > 0) = gm.cluster_type(idx(idx > 0)); labels = types;
        for t = {'E', 'I', 'nonsomatic', 'noise'}
            fprintf('%s: %d units (%d good, %d mua)\n', t{1}, nnz(strcmp(types, t{1})), nnz(strcmp(types, t{1}) & strcmp(unit_labels, 'good')), ...
                nnz(strcmp(types, t{1}) & strcmp(unit_labels, 'mua')))
        end
    end
    if exist('plot_on', 'var') && plot_on
        plotWaveformClusters(idx, pcs, group, gm, basis, tWF, unit_labels, P, all_dir, types, stats)
    end
    return
end

n_clusters = gm.NumComponents;
[stats, goodIDs, mxWF, tWF, wf_shift, unit_labels] = GMM_stats(all_dir, option_only_good, stat_set, basis);
[idx, nlogL, P, logpdf, d2] = gm.cluster(table2array(stats));
warning('on',id)
% P: posterior probabilities of each Gaussian mixture component in gm given each observation in X
% logpdf: logarithm of the estimated probability density function (pdf) evaluated at each observation in X
% d2: squared Mahalanobis distance of each observation in X to each Gaussian mixture component in gm

% Categorize cells where 1. max is greater than min or 2. stats could not be
% calculated as unknown
mask1 =  (max(mxWF,[],2)>abs(min(mxWF,[],2))) | isnan(idx);
idx(mask1) = 0;

% Categorize cells where 3. trough is not at expected time (>2 sigma) as unknown
[~, mins] = min(mxWF,[],2);
tmins = tWF(mins);
if any(~isnan(wf_shift)) % aligned waveforms: undo the shift to get the original trough time
    shift = wf_shift(:)'; shift(isnan(shift)) = 0;
    tmins = tmins + shift*mean(diff(tWF));
end
mean_t = mean(tmins(idx>0));
std_t = std(tmins(idx>0));
% mask_tmin = abs(tmins)>0.25*1e-3; % (s)
mask_tmin = abs((tmins - mean_t)/std_t)>3;
idx(mask_tmin) = 0;

% Apply labels
labels = cell(size(idx));
labels(idx==0) = {'unknown'};
labels(idx==1) = {'E'};
labels(idx==2) = {'I'};
inds_valid = find(idx>0);

% Extract relevant post. prob and Mahalanobis distances, take log to plot
P_match = NaN(size(idx));
d2_match = NaN(size(idx));

P_match(inds_valid) = P(sub2ind(size(P),inds_valid, idx(inds_valid)));
P_plot = log(P_match);
d2_match(inds_valid) = d2(sub2ind(size(d2),inds_valid, idx(inds_valid)));
d2_plot = log(d2_match);
if n_clusters==2
  % Apply posterior probability criteria
    min_P = 0.7;
    mask = P_match<min_P & idx>0;
    idx(mask) = 3;
    labels(mask) = {'intermediate'};
    
    % Apply distance criteria
    max_d2a = 10;
    mask = d2_match>max_d2a & idx>0;
    idx(mask) = 4;
    labels(mask) = {'moderate_outlier'};
    
    
    max_d2b = 20;
    mask = d2_match>max_d2b & idx>0;
    idx(mask) = 5;
    labels(mask) = {'extreme_outlier'};
  
end
% function [idx,nlogL,P,logpdf,d2] = applyGMM(gm, stats, mxWF, tWF, plot_on)
fprintf('E %i, I %i, intermediate %i, moderate_outlier %i, extreme_outlier %i, unclassified %i\n', ...
    nnz(idx==1), nnz(idx==2), nnz(idx==3), nnz(idx==4), nnz(idx==5), nnz(idx==0))
% disp(mean(stats{idx==1,:}))
% disp(mean(stats{idx==2,:}))

if exist('plot_on','var') && plot_on
    %% Plot results
    fs = 1/mean(diff(tWF));
    fig = figure; % everything below draws into fig explicitly, even if another figure becomes current
    colormap(fig, cool)
    ax_names = {'', 'Ratio peak:trough deriv.', 'log squared Mahalanobis distance', 'log posterior P'};
    color_by = {idx, stats.pk_deriv, d2_plot, P_plot};
    hs = gobjects(1, 4);
    for k = [2 1 3 4] % categories on the left
        hs(k) = subplot(2, 4, k, 'Parent', fig);
        scatter3(hs(k), stats.duration/fs*1e3, stats.asymmetry, stats.rate_log10, 15, color_by{k}, 'filled'); axis(hs(k), 'vis3d')
        xlabel(hs(k), 'Duration tr-pk(ms)'); ylabel(hs(k), 'Asymmetry'); zlabel(hs(k), 'log10(rate)')
        c = colorbar(hs(k), 'north', 'AxisLocation', 'out'); c.Label.String = ax_names{k};
    end
    % discrete colors matching the waveform panels: unknown, E, I, intermediate, moderate, extreme outlier
    colormap(hs(1), [0.6 0.6 0.6; 1 0 0; 0 0 1; 0.6 0.4 1; 0.8 0.2 1; 0 0 0]); clim(hs(1), [-0.5 5.5])
    c = hs(1).Colorbar; c.Ticks = 0:5; c.TickLabelInterpreter = 'none';
    c.TickLabels = {'unknown','E','I','intermediate','moderate_outlier','extreme_outlier'};
    Link = linkprop(hs,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
    setappdata(fig, 'StoreTheLink', Link);
    axis(hs, 'square')
    set(hs, 'XLim', [0 1], 'YLim', [-1 1], 'ZLim', [-2.3 1.6], 'XTick', 0:.2:1)
    for h = hs; h.Position(4) = h.Position(4) - 0.05; end % lower the top edge (bottom unchanged) to leave room for the title

    % Plot waveforms
    for k = 1:4; hs(k) = subplot(2, 4, 4 + k, 'Parent', fig); hold(hs(k), 'on'); end
    if any(idx==1); plot(hs(1), tWF, zscore(mxWF(idx==1, :)'), 'k'); plot(hs(1), tWF, mean(zscore(mxWF(idx==1,:)'),2), 'r', 'LineWidth', 2); end
    title(hs(1), sprintf('Excitatory (n = %d)', sum(idx==1))); ylabel(hs(1), 'Z-scored amplitude')
    if any(idx==2); plot(hs(2), tWF, zscore(mxWF(idx==2, :)'), 'k'); plot(hs(2), tWF, mean(zscore(mxWF(idx==2,:)'),2), 'b', 'LineWidth', 2); end
    title(hs(2), sprintf('Inhibitory (n = %d)', sum(idx==2)))
    if any(idx==0); plot(hs(3), tWF, zscore(mxWF(idx==0, :)'), 'Color', [0.6 0.6 0.6]); end
    if any(idx==3); plot(hs(3), tWF, zscore(mxWF(idx==3, :)'), 'Color', [0.6 0.4 1]); end
    title(hs(3), {sprintf('\\color[rgb]{0.6 0.4 1}intermediate (n = %d)', sum(idx==3)), ...
        sprintf('\\color[rgb]{0.6 0.6 0.6}unknown (n = %d)', sum(idx==0))})
    if any(idx==4); plot(hs(4), tWF, zscore(mxWF(idx==4, :)'), 'Color', [0.8 0.2 1]); end
    if any(idx==5); plot(hs(4), tWF, zscore(mxWF(idx==5, :)'), 'Color', 'k'); end
    title(hs(4), {sprintf('\\color[rgb]{0.8 0.2 1}moderate outlier (n = %d)', sum(idx==4)), ...
        sprintf('\\color{black}extreme outlier (n = %d)', sum(idx==5))})
    axis(hs, 'square'); linkaxes(hs); set(hs, 'XLim', tWF([1 end]), 'YLim', [-6 3])
    grid(hs,'on')
    if isempty(option_only_good) || ~option_only_good
        unit_set = 'all units, good + MUA';
        if any(strcmp(unit_labels, 'noise')); unit_set = [unit_set ' + noise']; end
    elseif option_only_good == 2; unit_set = 'Kilosort''s original good units (training set)';
    else; unit_set = 'good units'; end
    [parents, sessions] = fileparts(cellstr(all_dir));
    is_sort = startsWith(sessions, 'kilosort'); % sort folder passed in: name it session + sort
    [~, parents] = fileparts(parents);
    sessions(is_sort) = strcat(cellstr(parents(is_sort)), {' '}, sessions(is_sort));
    annotation(fig, 'textbox', [0 0.93 1 0.04], 'String', sprintf('%s: GMM labels, %s (n = %d)', strjoin(cellstr(sessions), ', '), unit_set, numel(idx)), ...
        'Interpreter', 'none', 'FontSize', 12, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'EdgeColor', 'none')
    drawnow
end
end

function D = componentDensity(F, gm, nu)
% Density of each row of F under each component of gm: Gaussian, or multivariate t with nu degrees of
% freedom (location gm.mu, scale gm.Sigma) when nu is given
K = gm.NumComponents; d = size(F, 2); D = zeros(size(F, 1), K);
for k = 1:K
    S = gm.Sigma(:, :, min(k, end));
    if isempty(nu)
        D(:, k) = mvnpdf(F, gm.mu(k, :), S);
    else
        Xc = F - gm.mu(k, :); delta = sum((Xc / S) .* Xc, 2);
        D(:, k) = exp(gammaln((nu + d)/2) - gammaln(nu/2) - d/2*log(nu*pi) - 0.5*log(det(S)) - (nu + d)/2*log1p(delta/nu));
    end
end
end

function plotWaveformClusters(idx, pcs, group, gm, basis, tWF, unit_labels, P, all_dir, types, stats)
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
pc_name = '';
if type_view
    col_name = {'E', 'I', 'nonsomatic', 'noise'};
    col_mask = cellfun(@(c) strcmp(types(:), c), col_name, 'UniformOutput', false);
    colors = [0.85 0.1 0.1; 0.1 0.3 0.85; 0.1 0.6 0.2; 0.3 0.3 0.3];
    col_title = col_name;
    % Trough-to-peak width and asymmetry (as in the trough GMM) from the normalized waveforms
    w = basis.w_norm; dt_ms = mean(diff(tWF)) * 1e3; n = size(w, 1);
    [~, i_tr] = min(w, [], 2); width = NaN(n, 1); asym = NaN(n, 1);
    for i = 1:n
        [pk, i_pk] = max(w(i, i_tr(i):end)); pk_pre = max(w(i, 1:i_tr(i)));
        width(i) = (i_pk - 1) * dt_ms; asym(i) = (pk - pk_pre) / (pk + pk_pre);
    end
    % Noise by group (noise cluster + outliers): peak group black, trough group grey
    is_noise = strcmp(types(:), 'noise'); is_outlier = is_noise & group(:) == 1;   % "outlier" = drawn grey
    i_rate = size(basis.wf_coeff, 2) + 1;                                  % rate_log10's column in the z-scoring
    rate = stats.rate_log10 * basis.feat_sd(i_rate) + basis.feat_mu(i_rate);
    m1 = gm.group(1); S1 = ((w - m1.mu) * m1.coeff - m1.s_mu) ./ m1.s_sd;  % every unit on the trough-group PCs
    model_xyz = [stats.wf_pc1, S1(:, [1 3])];
else
    n_clusters = max(idx); col_mask = arrayfun(@(k) idx == k, 1:n_clusters, 'UniformOutput', false);
    colors = lines(n_clusters);
    k_gauss = arrayfun(@(m) m.gm.NumComponents, gm.group);                  % background (if any) follows its group's Gaussians
    has_bg = arrayfun(@(m) isfield(m, 'p0') && m.p0 > 0, gm.group);
    has_out = arrayfun(@(m) isfield(m, 'nu') && ~isempty(m.nu) && isfield(m, 'outlier_p') && ~isempty(m.outlier_p) && m.outlier_p > 0, gm.group);
    extra = has_bg | has_out; col_title = compose('C%d', 1:n_clusters);
    names = {' (background)', ' (outliers)'};
    if extra(1); col_title{k_gauss(1) + 1} = [col_title{k_gauss(1) + 1} names{1 + has_out(1)}]; end
    if extra(2); k2 = k_gauss(1) + extra(1) + k_gauss(2) + 1; col_title{k2} = [col_title{k2} names{1 + has_out(2)}]; end
end
n_col = numel(col_mask);
t_ms = tWF * 1e3;
lag_ms = sqrt(basis.acg_edges(1:end-1) .* basis.acg_edges(2:end)) * 1e3;   % geometric bin centers
fig = figure('Name', 'Whole-waveform GMM clusters');
% Fixed axes positions (not tiledlayout, which re-fits every panel when a 3D panel is rotated)
col_w = 0.95 / n_col; row_y = [0.68 0.43 0.05]; row_h = [0.2 0.16 0.3];  % room for titles and x labels
pos = @(row, k) [0.045 + (k - 1)*col_w, row_y(row), 0.78*col_w, row_h(row)];
pc_ax = gobjects(0); pc_group = [];
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
                sub = sprintf('%d peak (black), %d trough (grey)', nnz(in_k & ~is_outlier), nnz(in_k & is_outlier));
            else
                sub = sprintf('%d good (black), %d mua (grey), P<0.7: %d', nnz(in_k & is_good), nnz(in_k & ~is_good), ...
                    nnz(in_k & max(P, [], 2) < 0.7));
            end
            title(ax, {sprintf('%s: n = %d', col_title{k}, nnz(in_k)), sub}, 'FontWeight', 'normal')
            ylim(ax, [-1.2 1.2]); xlabel(ax, 'time (ms)'); if k == 1; ylabel(ax, 'waveform / |largest extremum|'); end
        else
            set(ax, 'XScale', 'log'); xlabel(ax, 'lag (ms)')
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
        scatter3(ax, pcs(rest, 1), pcs(rest, 2), pcs(rest, 3), 5, [0.7 0.7 0.7])
        scatter3(ax, pcs(mine & ~is_good, 1), pcs(mine & ~is_good, 2), pcs(mine & ~is_good, 3), 8, c_light)
        scatter3(ax, pcs(mine & is_good, 1), pcs(mine & is_good, 2), pcs(mine & is_good, 3), 8, c)
        if isempty(pc_name); pc_title = sprintf('group %d PCs 1, 3', g); else; pc_title = pc_name; end   % group view: x, y = group PC1, PC3; z = z feature
        title(ax, pc_title, 'FontWeight', 'normal')
        xlabel(ax, 'PC1'); ylabel(ax, 'PC3'); zlabel(ax, 'z feature')
        axis(ax, 'equal'); axis(ax, 'vis3d'); grid(ax, 'on')
        set(ax, 'XTickLabel', [], 'YTickLabel', [], 'ZTickLabel', [])           % no labels, so the panel can be bigger
        pc_ax(end + 1) = ax; pc_group(end + 1) = g; %#ok<AGROW>
    end
end
if type_view
    is_ei = ismember(types(:), {'E', 'I'}); is_cell = ~strcmp(types(:), 'noise');
    feat_xyz = [width, asym, rate];
    feat_lab = {'width (ms)', 'asymmetry', 'log_{10} rate (Hz)'}; model_lab = {'wf PC1 (split)', 'trough PC1', 'trough PC3'};
    panels = {feat_xyz, is_ei & is_good, false, 'good E/I', feat_lab; ...
              feat_xyz, is_ei, false, 'good+MUA E/I', feat_lab; ...
              model_xyz, is_cell, false, 'E/I/NS, model axes', model_lab; ...
              model_xyz, is_cell, true, '+ noise: peak (black), trough (grey)', model_lab};
    for k = 1:4
        [xyz, show, with_noise, ttl, lab] = panels{k, :};
        ax = axes(fig, 'Position', pos(3, k), 'PositionConstraint', 'innerposition'); view(ax, 3); hold(ax, 'on')
        xlabel(ax, lab{1}); ylabel(ax, lab{2}); zlabel(ax, lab{3})
        if with_noise   % noise by group: trough grey, peak black
            nz = is_noise & is_outlier; scatter3(ax, xyz(nz, 1), xyz(nz, 2), xyz(nz, 3), 5, [0.6 0.6 0.6])
            nz = is_noise & ~is_outlier; scatter3(ax, xyz(nz, 1), xyz(nz, 2), xyz(nz, 3), 5, [0 0 0])
        end
        for t = 1:3   % E, I, nonsomatic: MUA lighter, then good in the type color
            in_t = show & col_mask{t}; c = colors(t, :);
            scatter3(ax, xyz(in_t & ~is_good, 1), xyz(in_t & ~is_good, 2), xyz(in_t & ~is_good, 3), 8, c + 0.6*(1 - c))
            scatter3(ax, xyz(in_t & is_good, 1), xyz(in_t & is_good, 2), xyz(in_t & is_good, 3), 8, c)
        end
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
[parents, sessions] = fileparts(cellstr(all_dir));
[~, parents] = fileparts(parents);
n_pc = [gm.group.n_pc];
sgtitle(fig, sprintf('Two-step whole-waveform GMM, %d clusters (group 1: %d PCs, group 2: %d PCs): %s %s (n = %d)', max(idx), n_pc, ...
    strjoin(unique(parents), ', '), strjoin(unique(sessions), ', '), numel(idx)), 'Interpreter', 'none', 'FontSize', 10)
drawnow
end

