function applyQualityMetrics(ks_dir, min_type_p, min_spikes_good, min_fit_p)
%% Apply basic sorting quality metrics to uncurated sessions
% Make sure to record if manually sorted in the excel table!!!
% Cell type model (openEphys pipeline): good = Kilosort 4 called it good (its own contamination test:
%   ContamPct < 20% and a significantly empty refractory window), more than min_spikes_good spikes, cell
%   type E / I / nonsomatic, cell type posterior >= min_type_p, and for E / I units, fit to their E / I
%   component p >= min_fit_p (gmm_result.mat fit_p); noise = cell type noise; the rest mua.
%   cluster_ContamPct.tsv keeps Kilosort 4's values.
% Trough GMM labels (Intan pipeline): also contamination (wvStruct.contam, Kilosort 2's measure from
%   getSessionWaveforms) < 10% and no all-channel artifact (see below); cluster_ContamPct.tsv gets
%   wvStruct.contam.
% min_type_p (optional, cell type model only): omitted = no posterior rule
% min_spikes_good (optional, default 50)
% min_fit_p (optional, cell type model with an E / I split only): omitted or 0 = no fit rule

if ~exist('min_spikes_good', 'var') || isempty(min_spikes_good); min_spikes_good = 50; end
max_contam_good = 0.10; % trough GMM labels only
max_contam_mua = 1;   % trough GMM labels only (see below)
min_spikes_mua = 20;  % trough GMM labels only: less than N spikes = noise

% Find the KS directory
disp(ks_dir)

% Never relabel a phy-curated session (phy writes cluster_info.tsv when curation is saved)
if isfile(fullfile(ks_dir, 'cluster_info.tsv'))
    warning('%s has been curated in phy (cluster_info.tsv exists); not applying quality metrics', ks_dir)
    return
end

% Get the contam percents calculated in getSessionWaveforms
wvStruct = getfield(load(fullfile(ks_dir, 'waveformStruct.mat')), 'wvStruct');
gm_result = load(fullfile(ks_dir,'gmm_result.mat')); % labels: cell type of each unit (E / I / nonsomatic / noise), from step 6

% Delete old phy log files -- will interfere with update of these labels
if exist(fullfile(ks_dir,'phy.log'),'file')
    warning('Deleting phy logs? Enter any key to continue')
    keyboard
    rmdir(fullfile(ks_dir,'.phy'),'s')
    delete(fullfile(ks_dir,'phy.log'))

end

% Copy original files
if ~exist(fullfile(ks_dir, 'cluster_ContamPct_orig.tsv'),'file')
    copyfile(fullfile(ks_dir, 'cluster_ContamPct.tsv'), fullfile(ks_dir, 'cluster_ContamPct_orig.tsv'))
    copyfile(fullfile(ks_dir, 'cluster_KSLabel.tsv'), fullfile(ks_dir, 'cluster_KSLabel_orig.tsv'))
end

% Read the old cluster IDs and labels
T_cluster_ContamPct_old = readtable(fullfile(ks_dir, 'cluster_ContamPct_orig.tsv'),'FileType','Text');
cluster_id = T_cluster_ContamPct_old.cluster_id;
T_cluster_KSLabel_old = readtable(fullfile(ks_dir, 'cluster_KSLabel_orig.tsv'),'FileType','Text');
mask_good_orig = strcmp(T_cluster_KSLabel_old.KSLabel,'good');

if isfield(gm_result, 'cell_type_model')
    % Cell type model (openEphys pipeline, step 6) decides noise: units typed noise are noise; E, I and
    % nonsomatic units are good (Kilosort called them good, enough spikes, and the model assigned them
    % clearly: posterior of their type >= min_type_p, if given) or mua. No separate contamination check
    % (Kilosort 4's good label has one) and no artifact check (all-channel artifacts are typed noise).
    mask_good = mask_good_orig & wvStruct.nSpikes(:) > min_spikes_good & ismember(gm_result.labels(:), {'E', 'I', 'nonsomatic'});
    if exist('min_type_p', 'var') && ~isempty(min_type_p); mask_good = mask_good & gm_result.type_p(:) >= min_type_p; end
    if exist('min_fit_p', 'var') && ~isempty(min_fit_p) && min_fit_p > 0 && isfield(gm_result, 'fit_p')
        mask_good = mask_good & ~(gm_result.fit_p(:) < min_fit_p);   % NaN (not E / I, or no E / I split): no rule
    end
    mask_noise = strcmp(gm_result.labels(:),'noise');
else
    % Trough GMM labels (Intan pipeline): good also needs low contamination and no artifact; noise from
    % contamination, spike count, outliers and artifacts. Never moves a unit from mua to good.
    mask_good = wvStruct.contam(:)<max_contam_good & ~isnan(wvStruct.contam(:)) & wvStruct.nSpikes(:)>min_spikes_good;
    max_all = squeeze(max(wvStruct.waveFormsMean,[],1));     % artifact: waveform high on every channel
    max_range = (max(max_all) - min(max_all))./max(max_all);
    max_range_threshold = 0.5;
    mask_artifact = max_range(:)<max_range_threshold;
    mask_noise = wvStruct.contam(:)>max_contam_mua | wvStruct.nSpikes(:)<min_spikes_mua;
    mask_good = mask_good & ~mask_noise & ismember(gm_result.labels(:),{'E','I','intermediate','moderate_outlier','unknown'}) & ~mask_artifact & mask_good_orig;
    mask_noise = mask_noise | ismember(gm_result.labels(:),{'extreme_outlier'}) | mask_artifact;
end
mask_mua = ~mask_noise & ~mask_good;
if nnz(mask_mua)+nnz(mask_good)+nnz(mask_noise) ~= length(mask_noise); error('check masks'); end



% TODO: Label units according to firing rate stability

%{
    figure;  % PLOT results
    subplot(2,2,1); plot(wvStruct.mxWF(mask_good&strcmp(gm_result.labels,'E'),:)','r');  title({T.filename{ii}, 'E'},'Interp','none'); grid on
    subplot(2,2,2); plot(wvStruct.mxWF(mask_good&strcmp(gm_result.labels,'I'),:)','b');  title('I'); grid on
    subplot(2,2,3); plot(wvStruct.mxWF(mask_mua,:)','k');  title('MUA'); grid on
    subplot(2,2,4); plot(wvStruct.mxWF(mask_noise,:)','Color',.3*[1 1 1]); title('Noise'); grid on
    linkaxes; ylim([-1.2 .6]*1e3);
%}

% Save new "good", "mua", "noise" labels
KSLabel = cell(length(cluster_id),1);
KSLabel(mask_good) = deal({'good'});
KSLabel(mask_mua) = deal({'mua'});
KSLabel(mask_noise) = deal({'noise'});
fprintf('ngood new %i, ngood orig %i\n', nnz(mask_good), nnz(strcmp(T_cluster_KSLabel_old.KSLabel,'good')))
T_cluster_KSLabel_new = table(cluster_id, KSLabel);
group = KSLabel;
T_cluster_group_new = table(cluster_id, group);
writetable(T_cluster_KSLabel_new, fullfile(ks_dir, 'cluster_KSLabel.tsv'),'FileType','Text','Delimiter','tab','WriteVariableNames',true);
writetable(T_cluster_group_new, fullfile(ks_dir, 'cluster_group.tsv'),'FileType','Text','Delimiter','tab','WriteVariableNames',true);

% cluster_ContamPct.tsv: Kilosort 4's own values for the cell type model (restored, in case an earlier run
% overwrote them); otherwise wvStruct.contam, the contamination used above, so Phy sees the same numbers
if isfield(gm_result, 'cell_type_model')
    copyfile(fullfile(ks_dir, 'cluster_ContamPct_orig.tsv'), fullfile(ks_dir, 'cluster_ContamPct.tsv'))
else
    ContamPct = round(wvStruct.contam*100 *10)/10;
    T_cluster_ContamPct_new = table(cluster_id,ContamPct);
    writetable(T_cluster_ContamPct_new, fullfile(ks_dir, 'cluster_ContamPct.tsv'),'FileType','Text','Delimiter','tab','WriteVariableNames',true); % TODO change from temp
end

% Update wvStruct
wvStruct.goodLabels = KSLabel;
wvStruct.goodIDs = wvStruct.goodIDs(:);
save(fullfile(ks_dir, 'waveformStruct.mat'),'wvStruct')

