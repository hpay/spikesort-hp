function [gm, basis] = fitCellTypeModel(ks_subdir, opts)
% FITCELLTYPEMODEL  Fit the pooled cell type model (E / I / nonsomatic / noise) used by step 6 of
% runSingleKilosortOpenEphys.m, and plot it.
%
%   fitCellTypeModel                                  % defaults below: saves src\cell_type_model_VW.mat
%   fitCellTypeModel('kilosort2_output', 'model_file', fullfile(code_dir, 'results', 'model_ks2.mat'))
%
% Whole-waveform features from each sort's waveformStruct.mat (cellTypeFeatures 'waveform'): peak-channel mean
% waveform (aligned on its largest |extremum|, normalized, sign kept), log firing rate, spatial spread.
% 1. makeCellTypeModel: split units into trough- and peak-dominated groups on waveform PC1, then refit a waveform
%    PCA within each group (3 PCs).
% 2. Within each group, a mixture of multivariate t-distributions (nu = 4; 3 trough + 2 peak clusters) on
%    the group PCs + log rate + spread (features per group in other). Units whose fit to their own cluster
%    has p < 0.05 (squared distance / d ~ F(d, nu)) are outliers.
%    Plots each cluster, then each cell type.
% 3. Cell type of each cluster: outliers and clusters whose mean ACG over 1-1.5 ms (log10 ratio to Poisson
%    + 0.1) is above 0 (no refractory dip) = noise; the remaining trough clusters are somatic; remaining
%    peak clusters = nonsomatic.
% 4. E vs I among the somatic units: 2-component t-mixture on trough-to-peak width, asymmetry and log rate, fit on
%    good somatic units only and applied to good and mua; wider component = E (gm.ei, applied by
%    cellTypesFromClusters).
% Steps 1-3 use all units of the sorts (each needs waveformStruct.mat, i.e. step 5 of the pipeline), so they
% don't depend on any good / mua / noise labels; step 4 uses the current good labels (cluster_group.tsv,
% step 7). Saves gm and basis to model_file (default src\cell_type_model_<region>.mat; one model per brain
% region, VW = visual wulst).

arguments
    ks_subdir = 'kilosort4_mirrored_interleaved_nodrift'   % sort folder in each session
    opts.sessions = {'MH02_20261111', 'MH02_20261118', 'MH02_20261121'}   % sessions to train on
    opts.region = 'VW'              % brain region the sessions are from; names the default model file
    opts.data_root = 'D:\dataMH'
    opts.model_file = ''            % default: src\cell_type_model_<region>.mat, read by step 6 of runSingleKilosortOpenEphys.m
end
ks_dirs = fullfile(opts.data_root, opts.sessions, ks_subdir);
code_dir = fileparts(fileparts(mfilename('fullpath')));
model_file = opts.model_file;
if isempty(model_file); model_file = fullfile(code_dir, 'src', ['cell_type_model_' opts.region '.mat']); end
units = 0;                          % all units (see header)
n_clusters = [3 2];                 % trough, peak
nu = 4; outlier_p = 0.05;
other = {'rate_log10', 'spread'};   % features besides the group PCs. Tried and rejected for the trough group:
                                    % adding width + asymmetry, and replacing the PCs with peak heights and times
                                    % (both split the E cells and lost the trough-noise cluster)
n_pc = 3;

addpath(genpath(fullfile(code_dir, 'src')))

%% 1. Fit: split on waveform PC1, waveform PCs within each group, t-mixture within each group
% (per-unit ACGs and spread are cached in each sort folder, so only the first fit on a sort computes them)
[gm, ~, ~, basis] = makeCellTypeModel(ks_dirs, units, n_clusters, 'waveform', [], other, n_pc, false, false, ...
    struct('nu', nu, 'outlier_p', outlier_p));

%% 2. Apply once: each unit's cluster, then the cell type of each cluster; plot and save
[idx, ~, ~, ~, P, res] = applyCellTypeModel(gm, ks_dirs, units, 0, 'waveform', basis);
gm.cluster_type = assignCellTypes(idx, res.basis, gm);
fprintf('cluster types: %s\n', strjoin(compose('C%d %s', (1:numel(gm.cluster_type))', string(gm.cluster_type')), ', '))

%% 3. E vs I: 2-component t-mixture on trough-to-peak width, asymmetry and log rate, fit on the good units of
% the somatic (E + I) clusters, then applied to all their units, good and mua (cellTypesFromClusters).
% fitTMixture numbers components by their mean on the first feature (width), so component 2 (wider) is E.
somatic = ismember(idx, find(ismember(gm.cluster_type, {'E', 'I'})));
fit_set = somatic & strcmp(res.unit_labels(:), 'good');
ei_features = {'width_ms', 'asym', 'rate_log10'};   % per-unit values in res.basis
X_ei = cell2mat(cellfun(@(f) res.basis.(f)(fit_set), ei_features, 'UniformOutput', false));
g_ei = fitTMixture(X_ei, 2, nu);
gm.ei = struct('gm', g_ei, 'nu', nu, 'e_comp', 2, 'features', {ei_features});
[res.types, res.type_p, res.ei_fit_p] = cellTypesFromClusters(gm, idx, P, res.basis);
fprintf('E / I fit on %d good somatic units, applied to all %d: E %d, I %d\n', nnz(fit_set), nnz(somatic), ...
    nnz(strcmp(res.types, 'E')), nnz(strcmp(res.types, 'I')))
plotCellTypeModel(gm, res, 3)                                   % each cluster, then each cell type
save(model_file, 'gm', 'basis', 'ks_dirs')
fprintf('saved %s\n', model_file)
end


function cluster_type = assignCellTypes(idx, b, gm)
% Cell type of each cluster from the units assigned to it (see the header for the rules); b = per-unit
% curves and shape from cellTypeFeatures (acg, width_ms)
lag_ms = sqrt(b.acg_edges(1:end-1) .* b.acg_edges(2:end)) * 1e3; early = lag_ms >= 1 & lag_ms <= 1.5;
k_g = arrayfun(@(m) m.gm.NumComponents, gm.group); has_out = arrayfun(@(m) m.outlier_p > 0, gm.group);
n_k = sum(k_g + has_out); cluster_type = repmat({'noise'}, 1, n_k); width = NaN(1, n_k);
outlier_k = [k_g(1) + 1, n_k]; outlier_k = outlier_k(has_out);
for k = 1:n_k
    in_k = idx == k;
    width(k) = median(b.width_ms(in_k));
    if ~ismember(k, outlier_k) && mean(b.acg(in_k, early), 'all') <= 0
        if k <= k_g(1); cluster_type{k} = 'somatic'; else; cluster_type{k} = 'nonsomatic'; end
    end
end
somatic = find(strcmp(cluster_type, 'somatic'));
[~, order] = sort(width(somatic), 'descend');
cluster_type(somatic(order(1))) = {'E'}; cluster_type(somatic(order(2:end))) = {'I'};
end
