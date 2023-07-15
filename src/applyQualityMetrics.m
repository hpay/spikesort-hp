function applyQualityMetrics(ks_dir)
%% Apply basic sorting quality metrics to uncurated sessions
% Make sure to record if manually sorted in the excel table!!!

max_contam_good = 0.10;
max_contam_mua = 1;
min_spikes_good = 50; % Label units with less than N spikes as mua
min_spikes_mua = 20;  % Label units with less than N spikes as noise

    % Find the KS directory
    disp(ks_dir)
    
    % Get the contam percents calculated in getSessionWaveforms
    wvStruct = getfield(load(fullfile(ks_dir, 'waveformStruct.mat')), 'wvStruct');
    gm_result = load(fullfile(ks_dir,'gmm_result.mat'));
    
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
    
    % Label units according to contamination
    mask_noise = wvStruct.contam>max_contam_mua | wvStruct.nSpikes<min_spikes_mua;
    mask_good = ~mask_noise & (wvStruct.contam<max_contam_good & ~isnan(wvStruct.contam) & wvStruct.nSpikes>min_spikes_good);
    
        % NEW: check if a waveform is high on every channel
    max_all = squeeze(max(wvStruct.waveFormsMean,[],1));
    max_range = (max(max_all) - min(max_all))./max(max_all);
    max_range_threshold = 0.5;
    mask_artifact = max_range(:)<max_range_threshold;
    
    % Label units according to GMM
    mask_good = mask_good & ismember(gm_result.labels,{'E','I'}) & ~mask_artifact;
    mask_noise = mask_noise | ismember(gm_result.labels,{'extreme_outlier','unknown'}) | mask_artifact;
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
    writetable(T_cluster_KSLabel_new, fullfile(ks_dir, 'cluster_KSLabel.tsv'),'FileType','Text','Delimiter','tab','WriteVariableNames',true);
    writetable(T_cluster_KSLabel_new, fullfile(ks_dir, 'cluster_group.tsv'),'FileType','Text','Delimiter','tab','WriteVariableNames',true);
    
    % Update cluster_ContamPct.tsv so Phy sees the same numbers
    ContamPct = round(wvStruct.contam*100 *10)/10;
    T_cluster_ContamPct_new = table(cluster_id,ContamPct);
    writetable(T_cluster_ContamPct_new, fullfile(ks_dir, 'cluster_ContamPct.tsv'),'FileType','Text','Delimiter','tab','WriteVariableNames',true); % TODO change from temp
    
    % Update wvStruct
    wvStruct.goodLabels = KSLabel;
    wvStruct.goodIDs = wvStruct.goodIDs(:);
    save(fullfile(ks_dir, 'waveformStruct.mat'),'wvStruct')

