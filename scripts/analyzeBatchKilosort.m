function analyzeBatchKilosort(T, data_dir, ksDirFun, results_dir, overwrite_gmm)


%% Get waveforms for all sessions in table T
only_good = 0; % Process all for now, select good later
for ii = height(T):-1:1
    
    % Find the binary directory
    root_dir = fullfile(data_dir, T.filename{ii});
    raw_dir_temp = dir(fullfile(root_dir, 'raw*'));
    if isempty(raw_dir_temp); continue; end
    raw_dir = fullfile(raw_dir_temp.folder, raw_dir_temp.name);
    ks_dir = ksDirFun(root_dir);

    % Check if already done with the latest sorting unit labels
    if exist(fullfile(ks_dir, 'waveformStruct.mat'),'file')
        wvStruct = getfield(load(fullfile(ks_dir, 'waveformStruct.mat')),'wvStruct');
        
        % Get the latest phy labels
        [unit_ID,cluster_labels] = getPhyClusterLabels(ks_dir);
        
        % Check if everything is already identical
        if length(unit_ID)==length(wvStruct.goodIDs) && all(unit_ID(:)==wvStruct.goodIDs(:)) && all(strcmp(wvStruct.goodLabels,cluster_labels)) && isfield(wvStruct, 'nSpikes')
            continue;
            
            % Check if units are identical, but labels (good/mua/noise) need to be updated
        elseif length(unit_ID)==length(wvStruct.goodIDs) && all(unit_ID(:)==wvStruct.goodIDs(:)) && ~all(strcmp(wvStruct.goodLabels,cluster_labels)) 
            % Re-save in case any updates to KS labels
            wvStruct.goodLabels = cluster_labels;
            save(fullfile(ks_dir,'waveformStruct.mat'),'wvStruct')
            fprintf('Saved new waveform labels %s\n', ks_dir)
            continue;
        end
        
    end
    
    % Save waveforms for cell type clustering
    disp(T(ii,:))
    wvStruct = getSessionWaveforms(raw_dir, ks_dir, only_good);
    save(fullfile(ks_dir, 'waveformStruct.mat'), 'wvStruct','-v7.3')
    fprintf('Saved new waveforms %s\n', ks_dir)

end

%% Generate GMM based on all curated sessions - rerun after sorting new sessions
gmm_name = 'latestGMM.mat';
if ~exist(fullfile(results_dir, gmm_name),'file') || overwrite_gmm
    T_load = T(strcmp(T.manually_sorted,'yes'),:);
    option_only_good = 1;
    n_clusters = 2;
    gm = GMM_make(cellfun(@(x) fullfile(data_dir, x), T_load.filename,'Uni',0), option_only_good, n_clusters);
    save(fullfile(results_dir, gmm_name),'gm','T')
end


%% Apply GMM results to all sessions
gm = getfield(load(fullfile(results_dir, gmm_name)),'gm');
option_only_good = 0;
plot_on = 1;
for ii = height(T):-1:1
    ks_dir = ksDirFun(fullfile(data_dir, T.filename{ii}));
    if ~exist(ks_dir,'dir'); continue; end

    if ~exist(fullfile(ks_dir,'gmm_result.mat'),'file') || overwrite_gmm
        disp(T.filename{ii})
        [id,labels, stats] = GMM_apply(gm, fullfile(data_dir, T.filename{ii}), option_only_good, plot_on);
        save(fullfile(ks_dir,'gmm_result.mat'),'id','labels','stats')
        
        wvStruct = getfield(load(fullfile(ks_dir, 'waveformStruct.mat')), 'wvStruct');
        wvStruct.typeLabels = labels;
        save(fullfile(ks_dir, 'waveformStruct.mat'),'wvStruct')
    end
end

%% Apply basic sorting quality metrics to uncurated sessions
% Make sure to record if manually sorted in the excel table!!!
for ii = height(T):-1:1
    ks_dir = ksDirFun(fullfile(data_dir, T.filename{ii}));
    if ~exist(ks_dir,'dir'); continue; end

    if strcmp(T.manually_sorted{ii},'yes'); disp('Skipping, manually sorted');
        continue;
    end
    applyQualityMetrics(ks_dir); % updates wvStruct.goodLabels in waveformStruct.mat
end

