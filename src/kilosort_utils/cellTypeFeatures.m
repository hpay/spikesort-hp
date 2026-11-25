function [stats, goodIDs, mxWF, tWF, wf_shift, unit_labels, basis] = cellTypeFeatures(all_dir, option_only_good, stat_set, basis)
% wf_shift: alignment shift (samples) from getSessionWaveforms(..., align=true); NaN if not aligned
% unit_labels: each unit's current good/mua/noise label (wvStruct.goodLabels)
% option_only_good: 0/[] all units, 1 good, 2 Kilosort's original good, or a cell array of labels to
%   keep, e.g. {'good','mua'}
% stat_set: 'trough' (default) = rate_log10, duration, asymmetry, pk_deriv (trough-based, for E/I)
%   'waveform' = whole-waveform features that also work for positive-going units:
%       wf_pc1..6: PCs of the peak-channel mean waveform, aligned on its largest |extremum| and
%           divided by that |extremum| (sign kept). The basis is fit on the above-median-amplitude
%           units unless one is passed in (basis, as returned here).
%       rate_log10
%       acg_pc1..3: PCs of the autocorrelogram (log-spaced bins 1-500 ms; log10 of the ratio to a Poisson
%           process at the unit's rate, +0.1). A basis passed in keeps its own bins and offset (some
%           older models: 0.5-500 ms or 3-500 ms with +0.5).
%       spread: mean peak-to-peak amplitude of channels 25-60 um from the peak channel, as a
%           fraction of the peak channel's (spatial fall-off)
%       All 11 are z-scored (mean and SD stored in basis).
%       pk1_height, pk1_time, pk2_height, pk2_time: largest value of the normalized waveform before the
%           trough (peak 1) and after it (peak 2), and their times from the trough (ms); z-scored
%           separately (basis.peak_mu / peak_sd). basis.width_ms (= peak 2 time) and basis.asym
%           ((pk2 - pk1) / (pk2 + pk1)) keep the unscaled trough-to-peak width and asymmetry for plots.
% basis: waveform and ACG PCA bases ('waveform' only), to score other units the same way
if ~exist('stat_set', 'var') || isempty(stat_set); stat_set = 'trough'; end
is_wf = strcmp(stat_set, 'waveform');
acg_edges = logspace(log10(1e-3), log10(0.5), 31); % s
acg_offset = 0.1;
if exist('basis', 'var') && ~isempty(basis)       % score with the bins and offset the basis was fit with
    acg_edges = basis.acg_edges;
    if isfield(basis, 'acg_offset'); acg_offset = basis.acg_offset; else; acg_offset = 0.1; end
end
acg = []; spread = [];

if istable(all_dir)
    R = all_dir;
    meanRate = R.rate;
    mxWF = cell2mat([R.waveform]')';
    goodIDs = ones(size(meanRate));
    wf_shift = NaN(size(meanRate));
    unit_labels = repmat({''}, size(meanRate));
    dir_idx = ones(size(meanRate));
    fs = 1e6;
    tWF = ((0:size(mxWF,2)-1)/fs - 3.08e-4)*1e3;
else

    if ~iscell(all_dir); all_dir = {all_dir}; end

    meanRate = [];
    dir_idx = [];
    mxWF = []; % cat(1,theseWaveforms.mxWF);
    goodIDs = []; %cat(1,theseWaveforms.goodIDs);
    wf_shift = [];
    unit_labels = {};

    for ii = 1:length(all_dir)
        ks_dir = all_dir{ii}; % a sort folder (name contains kilosort2/kilosort4), or a session folder containing kilosort2_output
        [~, dir_name] = fileparts(ks_dir);
        if ~contains(dir_name, 'kilosort'); ks_dir = fullfile(ks_dir, 'kilosort2_output'); end
        wvStruct = getfield(load(fullfile(ks_dir,'waveformStruct.mat')),'wvStruct');

        % Only keep "good" clusters. If sorted with Phy, this will be the final sorting
        % option_only_good = 2: Kilosort's own "good" labels, before applyQualityMetrics/phy
        if iscell(option_only_good)
            mask = ismember(wvStruct.goodLabels(:), option_only_good);
        elseif ~isempty(option_only_good) && option_only_good == 2
            ks_label_file = fullfile(ks_dir, 'cluster_KSLabel_orig.tsv'); % saved by applyQualityMetrics
            if ~isfile(ks_label_file); ks_label_file = fullfile(ks_dir, 'cluster_KSLabel.tsv'); end
            T_ks = readtable(ks_label_file, 'FileType', 'text');
            mask = ismember(wvStruct.goodIDs(:), T_ks.cluster_id(strcmp(T_ks.KSLabel, 'good')));
        elseif ~isempty(option_only_good) && option_only_good
            mask = strcmp(wvStruct.goodLabels,'good');
        else
            mask = true(size(wvStruct.max_site));
        end
        if is_wf; mask = mask(:) & ~any(isnan(wvStruct.mxWF), 2); end % units with no waveform

        meanRate = cat(1, meanRate, wvStruct.meanRate(mask,:));
        mxWF = cat(1,mxWF, wvStruct.mxWF(mask,:));
        goodIDs = cat(1,goodIDs, wvStruct.goodIDs(mask(:)));
        if isfield(wvStruct, 'wf_shift'); shift_i = wvStruct.wf_shift(mask(:)); else; shift_i = NaN(nnz(mask),1); end
        wf_shift = cat(1, wf_shift, shift_i);
        unit_labels = cat(1, unit_labels, wvStruct.goodLabels(mask(:)));
        dir_idx = cat(1, dir_idx, repmat(ii, nnz(mask), 1));   % which sort folder each unit is from

        if is_wf % autocorrelograms and spatial spread, cached in the sort folder
            feat = unitFeatures(ks_dir, wvStruct, acg_edges);
            acg = cat(1, acg, log10(feat.acg_ratio(mask(:), :) + acg_offset));
            spread = cat(1, spread, feat.spread(mask(:)));
        end
    end

    fs = size(wvStruct.mxWF,2)/wvStruct.spkDur;

    tWF = linspace(-wvStruct(1).spkOffset,-wvStruct(1).spkOffset+wvStruct(1).spkDur,...
        fs*wvStruct(1).spkDur);

end

if is_wf
    %% Whole-waveform stats
    w = mxWF - median(mxWF(:, 1:10), 2);              % baseline from the first 0.33 ms
    amp = max(abs(w), [], 2);
    w_norm = w ./ amp;                                % largest extremum = +1 or -1
    new_basis = ~exist('basis', 'var') || isempty(basis);
    if new_basis
        fit_set = amp > median(amp);
        [basis.wf_coeff, ~, ~, ~, basis.wf_explained, basis.wf_mu] = pca(w_norm(fit_set, :), 'NumComponents', 6);
        [basis.acg_coeff, ~, ~, ~, basis.acg_explained, basis.acg_mu] = pca(acg, 'NumComponents', 3);
    end
    wf_pcs = (w_norm - basis.wf_mu) * basis.wf_coeff;
    acg_pcs = (acg - basis.acg_mu) * basis.acg_coeff;
    X = [wf_pcs, log10(meanRate), acg_pcs, spread];
    if new_basis; basis.feat_mu = mean(X, 'omitnan'); basis.feat_sd = std(X, 'omitnan'); end
    X = (X - basis.feat_mu) ./ basis.feat_sd;          % z-scored, so each feature starts on the same scale
    stats = array2table(X, 'VariableNames', [compose('wf_pc%d', 1:size(wf_pcs, 2)), {'rate_log10'}, ...
        compose('acg_pc%d', 1:size(acg_pcs, 2)), {'spread'}]);
    % Trough-to-peak width (ms) and asymmetry of the normalized waveform, z-scored on their own so models
    % fit before these existed keep their feature scaling
    peaks = troughPeaks(w_norm, mean(diff(tWF)) * 1e3);   % pk1 height, pk1 time, pk2 height, pk2 time
    if ~isfield(basis, 'peak_mu'); basis.peak_mu = mean(peaks); basis.peak_sd = std(peaks); end
    peak_z = (peaks - basis.peak_mu) ./ basis.peak_sd;
    stats.pk1_height = peak_z(:, 1); stats.pk1_time = peak_z(:, 2);
    stats.pk2_height = peak_z(:, 3); stats.pk2_time = peak_z(:, 4);
    width_ms = peaks(:, 4); asym = (peaks(:, 3) - peaks(:, 1)) ./ (peaks(:, 3) + peaks(:, 1));   % for plotting
    basis.acg_edges = acg_edges; basis.acg_offset = acg_offset;
    basis.acg = acg; basis.w_norm = w_norm;           % per-unit curves and shape, for plotting
    basis.width_ms = width_ms; basis.asym = asym;
    basis.rate_log10 = log10(meanRate);               % unscaled log10 rate (Hz)
    basis.dir_idx = dir_idx;                          % index into all_dir of each unit's sort folder
    return
end
basis = [];

%% Compute stats
mxWF_d = diff(mxWF')';
numUnits = length(meanRate);
stats = array2table(NaN(numUnits,4)); stats.Properties.VariableNames = {'rate_log10','duration','asymmetry','pk_deriv'};
for nUnit = 1:numUnits
    % get trough
    [tr, trInd] = min(mxWF(nUnit,:));

    % get peak (after trough)
    [pk, pkInd] = max(mxWF(nUnit, trInd:end));
    pkInd = pkInd + trInd - 1;

    % get pre-peak (before trough)
    pk_pre = max(mxWF(nUnit, 1:trInd));

    % get derivative trough and peak
    [tr_d,trInd_d] = min(mxWF_d(nUnit,:));
    [pk_d,pkInd_d] = max(mxWF_d(nUnit,:));

    % concatenate the following statS
    stats.rate_log10(nUnit) = log10(meanRate(nUnit));
    stats.duration(nUnit) =  pkInd-trInd;               % duration
    stats.asymmetry(nUnit) = (pk-pk_pre)./(pk+pk_pre);  % asymmetry
    stats.pk_deriv(nUnit) = log10(abs(pk_d./tr_d));     % pk derivative ratio
end


end

function feat = unitFeatures(ks_dir, wvStruct, acg_edges)
% Autocorrelogram (ratio to a Poisson process at the unit's rate, in acg_edges bins) and spatial spread of
% every unit in wvStruct order (NaN for units without a waveform). Saved in the sort folder as
% cell_type_features.mat and reloaded while the waveforms (mxWF), the spike assignments (spike_clusters.npy,
% which phy rewrites) and the bins are unchanged.
cache_file = fullfile(ks_dir, 'cell_type_features.mat');
if isfile(cache_file)
    feat = load(cache_file);
    d_cache = dir(cache_file); d_clu = dir(fullfile(ks_dir, 'spike_clusters.npy'));
    if isequal(feat.acg_edges, acg_edges) && isequaln(feat.mxWF, wvStruct.mxWF) && d_cache.datenum > d_clu.datenum
        return
    end
end
fs_spk = size(wvStruct.mxWF, 2) / wvStruct.spkDur;
spike_times = double(readNPY(fullfile(ks_dir, 'spike_times.npy'))) / fs_spk;
spike_clusters = readNPY(fullfile(ks_dir, 'spike_clusters.npy'));
% waveFormsMean columns are the sorted raw channels Kilosort used (bad ones excluded)
[chan_idx, order] = sort(double(readNPY(fullfile(ks_dir, 'channel_map.npy'))) + 1);
pos = readNPY(fullfile(ks_dir, 'channel_positions.npy')); pos = pos(order, :);
n = numel(wvStruct.goodIDs);
feat = struct('ids', wvStruct.goodIDs(:), 'acg_ratio', NaN(n, numel(acg_edges) - 1), 'spread', NaN(n, 1), ...
    'acg_edges', acg_edges, 'mxWF', wvStruct.mxWF);
for u = find(~any(isnan(wvStruct.mxWF), 2))'
    t = sort(spike_times(spike_clusters == wvStruct.goodIDs(u)));
    feat.acg_ratio(u, :) = autocorrelogram(t, acg_edges) ./ (numel(t) * wvStruct.meanRate(u) * diff(acg_edges));
    ptp = max(wvStruct.waveFormsMean(:, :, u)) - min(wvStruct.waveFormsMean(:, :, u));
    i_peak = find(chan_idx == wvStruct.max_site(u));
    dist = vecnorm(pos - pos(i_peak, :), 2, 2);
    feat.spread(u) = mean(ptp(dist >= 25 & dist <= 60)) / ptp(i_peak);
end
save(cache_file, '-struct', 'feat')
end

function peaks = troughPeaks(w, dt_ms)
% Each row of w: [peak 1 height, peak 1 time, peak 2 height, peak 2 time], peak 1 = largest value up to the
% trough (time <= 0), peak 2 = largest value from the trough on (time >= 0), times in ms from the trough
n = size(w, 1); [~, i_tr] = min(w, [], 2); peaks = NaN(n, 4);
for i = 1:n
    [pk1, i_pk1] = max(w(i, 1:i_tr(i))); [pk2, i_pk2] = max(w(i, i_tr(i):end));
    peaks(i, :) = [pk1, (i_pk1 - i_tr(i)) * dt_ms, pk2, (i_pk2 - 1) * dt_ms];
end
end

function counts = autocorrelogram(t, edges)
% Spike counts at positive lags (s) in each bin; t sorted
counts = zeros(1, numel(edges) - 1);
for k = 1:numel(t) - 1
    lag = t(1 + k:end) - t(1:end - k);
    if ~any(lag < edges(end)); break; end
    counts = counts + histcounts(lag, edges);
end
end
