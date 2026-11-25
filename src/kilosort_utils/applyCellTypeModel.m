function [idx, labels, stats, goodIDs, P, res] = applyCellTypeModel(gm, all_dir, option_only_good, plot_on, stat_set, basis)
% stat_set: 'trough' (default) or 'waveform' (see cellTypeFeatures); basis: from makeCellTypeModel ('waveform')
% 'waveform': gm is the two-step struct from makeCellTypeModel. idx = cluster number, numbered group 1 (negative
%   trough) first, then group 2 (0 if features missing); labels 'C1', 'C2', ..., or the cell type
%   (E / I / nonsomatic / noise; missing features = noise) if the model has a cluster_type map.
%   P: posterior of each cluster. res: everything plotCellTypeModel needs, so figures can be drawn
%   later without applying the model again, plus res.types and res.type_p (cellTypesFromClusters). plot_on: see plotCellTypeModel (1 per cell type, 2 per
%   cluster, 3 both)

if ~exist('option_only_good','var'); option_only_good = []; end
if ~exist('stat_set', 'var') || isempty(stat_set); stat_set = 'trough'; end
if ~exist('basis', 'var'); basis = []; end
id = 'stats:gmdistribution:cluster:MissingData'; warning('off',id)

if strcmp(stat_set, 'waveform')
    [stats, goodIDs, ~, tWF, ~, unit_labels, basis] = cellTypeFeatures(all_dir, option_only_good, stat_set, basis);
    group = 2 - (gm.split.cluster(stats.wf_pc1) == gm.split_low);
    has_bg = arrayfun(@(m) isfield(m, 'p0') && m.p0 > 0, gm.group);
    % t models: outlier_p > 0 moves units whose fit to their own component has p < outlier_p (squared
    % distance / d ~ F(d, nu)) to an outlier cluster, which takes the background's place
    has_out = arrayfun(@(m) isfield(m, 'nu') && ~isempty(m.nu) && isfield(m, 'outlier_p') && ~isempty(m.outlier_p) && m.outlier_p > 0, gm.group);
    k_group = arrayfun(@(m) m.gm.NumComponents, gm.group) + (has_bg | has_out);   % background / outliers: the group's last cluster
    idx = zeros(height(stats), 1); P = zeros(height(stats), sum(k_group)); pcs = NaN(height(stats), 3);
    for g = 1:2
        m = gm.group(g); in_g = group == g;
        pc_g = ((basis.w_norm(in_g, :) - m.mu) * m.coeff - m.s_mu) ./ m.s_sd;   % group PCs 1-3, z-scored
        other_g = gm.other; if isfield(m, 'other') && ~isempty(m.other); other_g = m.other; end   % per-group features, if set
        F = [pc_g(:, 1:m.n_pc), stats{in_g, other_g}];
        nu = []; if isfield(m, 'nu'); nu = m.nu; end                   % nu: t-distribution components
        if has_bg(g) || ~isempty(nu)   % posterior over the components (and the uniform background)
            L = (1 - m.p0) * m.gm.ComponentProportion .* mixtureDensity(F, m.gm, nu);
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
        z_label{g} = strrep(strrep(z_name, '_', ' '), 'rate log10', 'log_{10} rate');
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
    types = {}; type_p = []; ei_fit_p = [];
    if isfield(gm, 'cluster_type')   % model maps clusters to cell types: labels are E / I / nonsomatic / noise
        [types, type_p, ei_fit_p] = cellTypesFromClusters(gm, idx, P, basis); labels = types;
        for t = {'E', 'I', 'nonsomatic', 'noise'}
            fprintf('%s: %d units (%d good, %d mua)\n', t{1}, nnz(strcmp(types, t{1})), nnz(strcmp(types, t{1}) & strcmp(unit_labels, 'good')), ...
                nnz(strcmp(types, t{1}) & strcmp(unit_labels, 'mua')))
        end
    end
    res = struct('idx', idx, 'pcs', pcs, 'group', group, 'basis', basis, 'tWF', tWF, 'unit_labels', {unit_labels}, ...
        'all_dir', {all_dir}, 'stats', stats, 'z_label', {z_label}, 'ids', goodIDs, 'types', {types}, 'type_p', type_p, 'ei_fit_p', ei_fit_p);
    if exist('plot_on', 'var') && plot_on; plotCellTypeModel(gm, res, plot_on); end
    return
end

n_clusters = gm.NumComponents;
[stats, goodIDs, mxWF, tWF, wf_shift, unit_labels] = cellTypeFeatures(all_dir, option_only_good, stat_set, basis);
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
