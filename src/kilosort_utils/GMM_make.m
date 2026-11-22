function [gm, idx, labels, basis] = GMM_make(all_dir, option_only_good, n_clusters, stat_set, cov_type, other, n_pc_fixed, background)
% GMM_MAKE  Fit a Gaussian mixture to unit waveform features from one or more sort folders.
% 'trough' (default): E/I from trough features. Clusters are named by spike width: the widest mean
%   trough-to-peak duration is E, narrow is I. idx: 1 = E, 2 = I, 3 = noise (3 clusters).
% 'waveform': the two-step whole-waveform model (below); scripts/fitCellTypeModel.m builds the cell
%   type model on it. See https://www.mathworks.com/help/stats/gmdistribution.cluster.html for outputs.
% stat_set: 'trough' (default, E/I from trough features) or 'waveform' (whole-waveform features, see
%   GMM_stats; all units kept, number of clusters chosen by BIC if n_clusters is empty, clusters
%   numbered by their mean first waveform PC). basis: PCA bases to pass to GMM_apply ('waveform').
% cov_type ('waveform' only): 'full' (default), 'diagonal', or 'shared' (one full covariance shared by all
%   of a group's clusters, so no cluster can be much broader than the others); or a cell array with
%   one setting per group, e.g. {'shared', 'full'}
% other ('waveform' only): non-waveform features clustered with the group PCs, default
%   {'rate_log10', 'acg_pc1', 'acg_pc2', 'acg_pc3', 'spread'}
% n_pc_fixed ('waveform' only): number of group waveform PCs to cluster on; default [] = 2-3, enough for
%   90% of the group's variance
% background ('waveform' only, full covariance): true adds a uniform "catch-all" component to each group's
%   mixture (density = 1 / volume of the group's feature range, padded 10%), fit by EM together with the
%   Gaussians, so units that fit no cluster go to their own background cluster. Default false.
if ~exist('stat_set', 'var') || isempty(stat_set); stat_set = 'trough'; end
if ~exist('cov_type', 'var') || isempty(cov_type); cov_type = 'full'; end

%% Get waveforms from each session
[stats, goodIDs, mxWF, tWF, ~, ~, basis] = GMM_stats(all_dir, option_only_good, stat_set);

if strcmp(stat_set, 'waveform')
    % Two steps: 1. split units into two groups on wf_pc1 (negative vs positive main extremum) with a
    % 1-D, 2-component GMM. 2. Within each group, refit the waveform PCA on that group's normalized
    % waveforms, keep 2-3 PCs (enough for 90% of its variance), and fit a GMM on those PCs plus
    % the features in other (cluster count by BIC unless n_clusters = [n_group1 n_group2]).
    % gm is then a struct (split model + each group's basis and GMM) for GMM_apply.
    rng(2)
    if ~exist('other', 'var') || isempty(other); other = {'rate_log10', 'acg_pc1', 'acg_pc2', 'acg_pc3', 'spread'}; end
    X_other = stats{:, other};
    split = fitgmdist(stats.wf_pc1, 2, 'Replicates', 10);
    [~, split_low] = min(split.mu);
    group = 2 - (split.cluster(stats.wf_pc1) == split_low);  % 1 = low wf_pc1 (negative trough)
    gm = struct('split', split, 'split_low', split_low, 'other', {other});
    if ~iscell(cov_type); cov_type = {cov_type, cov_type}; end
    use_bic = ~(exist('n_clusters', 'var') && numel(n_clusters) == 2);
    if use_bic; fig = figure('Name', 'GMM BIC'); ax = axes(fig); hold(ax, 'on'); end
    for g = 1:2
        shared = strcmp(cov_type{g}, 'shared');
        fit_k = @(Xg, k) fitgmdist(Xg, k, 'CovarianceType', replace(cov_type{g}, 'shared', 'full'), ...
            'SharedCovariance', shared, 'Replicates', 20, 'RegularizationValue', 1e-3, 'Options', statset('MaxIter', 1000));
        [coeff, score, ~, ~, expl, mu_w] = pca(basis.w_norm(group == g, :));
        n_pc = min(3, max(2, find(cumsum(expl) >= 90, 1)));
        if exist('n_pc_fixed', 'var') && ~isempty(n_pc_fixed); n_pc = n_pc_fixed; end
        S = score(:, 1:3);                                    % 3 kept for plotting; the first n_pc clustered
        s_mu = mean(S); s_sd = std(S);
        Xg = [(S(:, 1:n_pc) - s_mu(1:n_pc)) ./ s_sd(1:n_pc), X_other(group == g, :)];
        if ~use_bic
            gm_g = fit_k(Xg, n_clusters(g));
        else
            k_try = 1:8; bic = NaN(size(k_try)); fits = cell(size(k_try));
            for k = k_try; fits{k} = fit_k(Xg, k); bic(k) = fits{k}.BIC; end
            [~, k_best] = min(bic); gm_g = fits{k_best};
            plot(ax, k_try, bic, 'o-', 'MarkerFaceColor', 'auto', 'DisplayName', ...
                sprintf('group %d (n = %d, %d PCs): minimum at %d', g, size(Xg, 1), n_pc, k_best))
        end
        [~, order] = sort(gm_g.mu(:, 1));                    % within a group, number clusters by mean group PC1
        if ~shared; Sigma = gm_g.Sigma(:, :, order); else; Sigma = gm_g.Sigma; end
        gm_g = gmdistribution(gm_g.mu(order, :), Sigma, gm_g.ComponentProportion(order));
        p0 = 0; bg_density = 0;
        if exist('background', 'var') && background
            bg_density = 1 / prod(1.2 * range(Xg));
            [gm_g, p0] = fitWithBackground(Xg, gm_g, bg_density, shared);
        end
        gm.group(g) = struct('coeff', coeff(:, 1:3), 'mu', mu_w, 's_mu', s_mu, 's_sd', s_sd, 'n_pc', n_pc, ...
            'explained', expl(1:3)', 'gm', gm_g, 'p0', p0, 'bg_density', bg_density);
    end
    if use_bic
        set(ax, 'FontSize', 10, 'TitleFontSizeMultiplier', 1, 'LabelFontSizeMultiplier', 1, 'Box', 'off', 'TickDir', 'out')
        xlabel(ax, 'number of clusters'); ylabel(ax, 'BIC'); legend(ax, 'Location', 'best', 'FontSize', 10)
        title(ax, sprintf('GMM within each wf PC1 group, group PCs (%s / %s covariance)', cov_type{:})); drawnow
    end
    [idx, labels] = GMM_apply(gm, all_dir, option_only_good, 1, stat_set, basis);
    return
end

% Train only on negative-trough units: trough features don't apply when the peak is bigger
% (GMM_apply labels those units unknown)
pos_peak = max(mxWF,[],2) > abs(min(mxWF,[],2));
fprintf('GMM_make: excluding %d of %d units with peak > trough from training\n', nnz(pos_peak), numel(pos_peak))
stats = stats(~pos_peak,:);

% Exclude very low rate cells from the template creation
% meanRateThresh = 0.05;  % min rate threshold
% mask = stats.rate_log10 < log10(meanRateThresh); 
% stats(mask,:) = [];
% mxWF(mask,:) = [];

if ~exist('n_clusters','var')
    n_clusters= 3;
end

%% Run GMM
rng(2)
X = table2array(stats);
gm = fitgmdist(X,n_clusters,'CovarianceType','diagonal','Replicates',25); % 
[idx,~,P,logpdf, d2] = gm.cluster(X);

% Identify clusters
mus = array2table(gm.mu); mus.Properties.VariableNames = stats.Properties.VariableNames;
[~, order] = sort(mus.duration, 'descend');  % E first = widest spike
sigmas = array2table(squeeze(gm.Sigma)'); sigmas.Properties.VariableNames = stats.Properties.VariableNames;
sigmas_test = table2array(sigmas(:,{'rate_log10','duration'}));

if n_clusters==3
    [~, ind_noise] = max(mean(sigmas_test./mean(sigmas_test),2));
    order(order==ind_noise)=[];    
    ind_exc = order(1);
    ind_inh  = find(~ismember(1:3, [ind_exc ind_noise]));
    ind_start = zeros(height(stats),1);
    ind_start(idx==ind_exc) = 1;
    ind_start(idx==ind_inh) = 2;
    ind_start(idx==ind_noise) = 3;
else
    ind_exc = order(1);    
    ind_start = zeros(height(stats),1);
    ind_start(idx==ind_exc) = 1;
    ind_start(idx~=ind_exc) = 2;
end

%% Re-run GMM with new initialization: 1 = E, 2 = I, 3 = unknown
% My previous notation: % 0 = I, 1 = E, NaN = unknown
cov_type = 'full'; % SC uses diagonal but I find clusters are not always round
gm = fitgmdist(X,n_clusters,'CovarianceType',cov_type,'Start',ind_start);

plot_on = 1;
[idx,labels] = GMM_apply(gm, all_dir, option_only_good, plot_on);

end

function [gm, p0] = fitWithBackground(X, gm, bg_density, shared)
% EM for a full-covariance Gaussian mixture plus a uniform background component of fixed density
% bg_density, starting from gm. shared: one covariance pooled across the Gaussians. Returns the refit
% Gaussians (proportions renormalized to sum to 1) and the background weight p0.
[n, d] = size(X); K = gm.NumComponents;
mu = gm.mu; Sigma = repmat(gm.Sigma, 1, 1, K / size(gm.Sigma, 3)); p = 0.95 * gm.ComponentProportion; p0 = 0.05;
loglik_old = -Inf;
for iter = 1:1000
    L = zeros(n, K + 1);
    for k = 1:K; L(:, k) = p(k) * mvnpdf(X, mu(k, :), Sigma(:, :, k)); end
    L(:, K + 1) = p0 * bg_density;
    loglik = sum(log(sum(L, 2)));
    R = L ./ sum(L, 2);                                   % responsibilities
    Nk = sum(R, 1);
    for k = 1:K
        mu(k, :) = R(:, k)' * X / Nk(k);
        Xc = X - mu(k, :);
        Sigma(:, :, k) = (Xc .* R(:, k))' * Xc / Nk(k) + 1e-3 * eye(d);  % same regularization as fitgmdist
    end
    if shared   % pool: weighted average of the clusters' covariances
        Sigma = repmat(sum(Sigma .* reshape(Nk(1:K), 1, 1, K), 3) / sum(Nk(1:K)), 1, 1, K);
    end
    p = Nk(1:K) / n; p0 = Nk(K + 1) / n;
    if loglik - loglik_old < 1e-6 * abs(loglik); break; end
    loglik_old = loglik;
end
if shared; Sigma = Sigma(:, :, 1); end
gm = gmdistribution(mu, Sigma, p / sum(p));
end

