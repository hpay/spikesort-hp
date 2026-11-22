%% Fit the pooled cell type model (E / I / nonsomatic / noise) used by step 6 of runSingleKilosortOpenEphys.m
% Whole-waveform features from each sort's waveformStruct.mat (GMM_stats 'waveform'): peak-channel mean
% waveform (aligned on its largest |extremum|, normalized, sign kept), log firing rate, spatial spread.
% 1. GMM_make: split units into trough- and peak-dominated groups on waveform PC1, then refit a waveform
%    PCA within each group (3 PCs). Also fits a Gaussian mixture, used here only for the split and bases.
% 2. Within each group, a mixture of multivariate t-distributions (nu = 4; 3 trough + 2 peak clusters) on
%    the group PCs + log rate + spread. Units whose fit to their own cluster has p < 0.05 (squared distance
%    / d ~ F(d, nu)) are outliers.
% 3. Cell type of each cluster: outliers and clusters whose mean ACG over 1-1.5 ms (log10 ratio to Poisson
%    + 0.1) is above 0 (no refractory dip) = noise; of the remaining trough clusters, longer median
%    trough-to-peak time = E, shorter = I; remaining peak clusters = nonsomatic.
% Fit on all units of the sorts below (each needs waveformStruct.mat, i.e. step 5 of the pipeline), so the
% fit doesn't depend on any good / mua / noise labels; saves gm and basis to src\cell_type_model_<region>.mat
% (one model per brain region; the committed cell_type_model_VW.mat is from visual wulst).
% The committed model was fit on the three MH02 sessions (good + MUA units, before step 7 used cell types);
% this script reproduces its cell type for 1085 of 1098 units.

%% Settings: sessions to train on (each needs waveformStruct.mat from step 5 of runSingleKilosortOpenEphys.m)
region = 'VW';                      % brain region the sessions are from (VW = visual wulst); names the model file
data_dir = 'D:\dataMH';
sessions = {'MH02_20261111', 'MH02_20261118', 'MH02_20261121'};
ks_subdir = 'kilosort4_mirrored_interleaved_nodrift';
ks_dirs = fullfile(data_dir, sessions, ks_subdir);
code_dir = fileparts(fileparts(mfilename('fullpath')));
model_file = fullfile(code_dir, 'src', ['cell_type_model_' region '.mat']);   % read by step 6 of runSingleKilosortOpenEphys.m
if exist('model_file_override', 'var'); model_file = model_file_override; end
units = 0;                          % all units (see header)
n_clusters = [3 2];                 % trough, peak
nu = 4; outlier_p = 0.05;

addpath(genpath(fullfile(code_dir, 'src')))

%% 1. Split and group bases (shared trough covariance + background, as in the model these settings came from)
[gm, ~, ~, basis] = GMM_make(ks_dirs, units, n_clusters, 'waveform', {'shared', 'full'}, {'rate_log10', 'spread'}, 3, true);

%% 2. t-mixture within each group
[stats, ~, ~, tWF, ~, ~, b] = GMM_stats(ks_dirs, units, 'waveform', basis);
group = 2 - (gm.split.cluster(stats.wf_pc1) == gm.split_low);
for g = 1:2
    m = gm.group(g); in_g = group == g;
    pc_g = ((b.w_norm(in_g, :) - m.mu) * m.coeff - m.s_mu) ./ m.s_sd;
    X = [pc_g(:, 1:m.n_pc), stats{in_g, gm.other}];
    gm.group(g).gm = fitTMixture(X, n_clusters(g), nu);
    gm.group(g).p0 = 0; gm.group(g).bg_density = 0; gm.group(g).nu = nu; gm.group(g).outlier_p = outlier_p;
end

%% 3. Cell type of each cluster, then plot and save
idx = GMM_apply(gm, ks_dirs, units, 0, 'waveform', basis);
gm.cluster_type = assignCellTypes(idx, b, basis, tWF, gm);
fprintf('cluster types: %s\n', strjoin(compose('C%d %s', (1:numel(gm.cluster_type))', string(gm.cluster_type')), ', '))
GMM_apply(gm, ks_dirs, units, 1, 'waveform', basis);
save(model_file, 'gm', 'basis', 'ks_dirs')
fprintf('saved %s\n', model_file)


function gm = fitTMixture(X, K, nu)
% EM for a mixture of K multivariate t-distributions with nu degrees of freedom (full covariances),
% starting from a full-covariance Gaussian fit; clusters numbered by their mean on the first feature
[n, d] = size(X);
rng(2); g0 = fitgmdist(X, K, 'CovarianceType', 'full', 'Replicates', 20, 'RegularizationValue', 1e-3, 'Options', statset('MaxIter', 1000));
[~, order] = sort(g0.mu(:, 1)); mu = g0.mu(order, :); Sigma = g0.Sigma(:, :, order); p = g0.ComponentProportion(order);
ll_old = -Inf;
for iter = 1:2000
    L = zeros(n, K); delta = zeros(n, K);
    for k = 1:K
        Xc = X - mu(k, :); delta(:, k) = sum((Xc / Sigma(:, :, k)) .* Xc, 2);
        L(:, k) = p(k) * exp(gammaln((nu + d)/2) - gammaln(nu/2) - d/2*log(nu*pi) - 0.5*log(det(Sigma(:, :, k))) ...
            - (nu + d)/2*log1p(delta(:, k)/nu));
    end
    ll = sum(log(sum(L, 2))); R = L ./ sum(L, 2);
    U = (nu + d) ./ (nu + delta);                            % downweights units far from a component
    for k = 1:K
        w = R(:, k) .* U(:, k);
        mu(k, :) = w' * X / sum(w);
        Xc = X - mu(k, :);
        Sigma(:, :, k) = (Xc .* w)' * Xc / sum(R(:, k)) + 1e-3 * eye(d);
    end
    p = mean(R, 1);
    if ll - ll_old < 1e-7 * abs(ll); break; end
    ll_old = ll;
end
fprintf('t-mixture (%d clusters): EM converged in %d iterations\n', K, iter)
gm = gmdistribution(mu, Sigma, p);
end

function cluster_type = assignCellTypes(idx, b, basis, tWF, gm)
% Cell type of each cluster from the units assigned to it (see the header for the rules)
lag_ms = sqrt(basis.acg_edges(1:end-1) .* basis.acg_edges(2:end)) * 1e3; early = lag_ms >= 1 & lag_ms <= 1.5;
dt_ms = mean(diff(tWF)) * 1e3;
k_g = arrayfun(@(m) m.gm.NumComponents, gm.group); has_out = arrayfun(@(m) m.outlier_p > 0, gm.group);
n_k = sum(k_g + has_out); cluster_type = repmat({'noise'}, 1, n_k); width = NaN(1, n_k);
outlier_k = [k_g(1) + 1, n_k]; outlier_k = outlier_k(has_out);
for k = 1:n_k
    in_k = idx == k; w = b.w_norm(in_k, :); [~, i_tr] = min(w, [], 2); tp = NaN(nnz(in_k), 1);
    for i = 1:nnz(in_k); [~, i_pk] = max(w(i, i_tr(i):end)); tp(i) = (i_pk - 1) * dt_ms; end
    width(k) = median(tp);
    if ~ismember(k, outlier_k) && mean(b.acg(in_k, early), 'all') <= 0
        if k <= k_g(1); cluster_type{k} = 'somatic'; else; cluster_type{k} = 'nonsomatic'; end
    end
end
somatic = find(strcmp(cluster_type, 'somatic'));
[~, order] = sort(width(somatic), 'descend');
cluster_type(somatic(order(1))) = {'E'}; cluster_type(somatic(order(2:end))) = {'I'};
end
