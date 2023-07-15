%% Script for batch processing kilosort spike sorting
% and backing up recently collected data to server!
%
% Instructions: Update RECORDING_DEPTH_CHICK.xlsx then run this script from
% within the scripts folder
%
% Necessary column names in spreadsheet:
%   filename: name of the session directory
%   bad_chan: list any bad channels in this format: [1, 6, 10]. If none: []
%   probe: name of probe channel map {H5, H6}
%   manually_sorted: 1 if manually curated in phy, 0 otherwise

T = readtable('C:\Users\Hannah\Dropbox\alab\Code\project2\data\RECORDING_DEPTH_CHICK.xlsx');
T = T(T.include>0,:);

% Re-spike sort and overwrite kilosort output?
overwrite_ks = 0;

% Re-calculate GMM based on manually curated sessions, and apply to all sessions?
overwrite_gmm = 0;

% Data folder
data_dir_local = 'D:\data';  % faster to process on locally saved files
data_dir_remote = 'Z:\Hannah\Ephys\Project2'; % Location to backup to

%% Set up

% Only process data that is actually present!
fnames = dir(data_dir_local);
fnames = {fnames.name};
T = T(ismember(T.filename, fnames),:);

% Save directory function -- ks results for each session saved here
ksDirFun = @(root_dir) fullfile(root_dir,'kilosort2_output');

% Repository folder
scripts_dir = fileparts(which(mfilename));
spikesort_hp_dir = fileparts(scripts_dir);
addpath(genpath(fullfile(spikesort_hp_dir,'src')))

% Add kilosort2.0 directory (from: https://github.com/MouseLand/Kilosort/releases/tag/v2.0)
code_dir = fileparts(spikesort_hp_dir); %
addpath(genpath(fullfile(code_dir, 'kilosort-2.0')))

% Common results folder (GMM etc)
results_dir = fullfile(fileparts(scripts_dir), 'results');
addpath(results_dir);

% warning("off","parallel:gpu:device:DeviceDeprecated");

%% Backup from local SSD to server
exceptions = {'kilosort2_output'}; % Back this up specifically later
runmode = 0;
for ii = 1:height(T)
    fprintf('\n%s\n',T.filename{ii})
    root_dir_local = fullfile(data_dir_local, T.filename{ii});  % e.g. D:\data\HC11_230129
    root_dir_remote = fullfile(data_dir_remote, T.filename{ii}); % e.g. Z:\Hannah\ephys\HC11_230129
    if ~exist(root_dir_local,'dir') || strcmp(root_dir_local, root_dir_remote); continue; end
    dirbackup(root_dir_local, root_dir_remote, runmode, exceptions)
end

%% Run Kilosort2
for ii = 1:height(T)
    % Root dir on server (Z:\Hannah\ephys\HC11_230129 etc)
    root_dir_local = fullfile(data_dir_local, T.filename{ii});
    
    % Find the binary directory
    raw_dir_temp = dir(fullfile(root_dir_local, 'raw*'));
    raw_dir = fullfile(raw_dir_temp.folder, raw_dir_temp.name);
    
    % Create save directory
    ks_dir = ksDirFun(root_dir_local);
    
    % If it exists already, either overwrite or skip
    if ~isempty(dir(fullfile(ks_dir,'*.npy')))
        if overwrite_ks
            fprintf('Overwriting %s\n', T.filename{ii})
            try
                rmdir(ks_dir,'s')
            catch
            end
        else
            fprintf('Already sorted %s\n', T.filename{ii})
            continue;
        end
    end
    mkdir(ks_dir)
    disp(T(ii,:))
    
    % set up options
    ops = [];
    ops.badChannels = eval(T.bad_chan{ii});
    ops.chanMap = [T.probe{ii} '.mat']; % Make sure config folder is on the search path
    ops = hp_config(ops);
    
    % Run and save results. Note: ops saved in rez.mat
    fprintf('\n Running Kilosort on directory %s \n', raw_dir),
    runSingleKilosort(raw_dir, ks_dir, ops);
    
    % Copy kilosort output to the server
    root_dir_remote = fullfile(data_dir_remote, T.filename{ii});
    ks_dir_remote = ksDirFun(root_dir_remote);
    copyfile(ks_dir, ks_dir_remote)
    
    % Edit params.py to reflect new amplifier.dat location
    fid = fopen(fullfile(ks_dir,'params.py'), 'r'); % Read params.py from local copy
    a = {}; while ~feof(fid); a{end+1} = fgets(fid);  pause(.0001);  end
    fclose(fid);
    a{1} = sprintf("dat_path = '%s'\n", strrep(fullfile(root_dir_remote, raw_dir_temp.name, 'amplifier.dat'),'\','/'));
    fid = fopen(fullfile(ks_dir_remote,'params.py'), 'w'); % Write correct params.py to remote server
    for jj = 1:length(a); fprintf(fid,a{jj});      end
    fclose(fid);
    fprintf('Results copied to server!\n')
    
end

%% Go manually label results in Phy!
% In Anaconda prompt:
%   conda activate phy_env
% Navigate to folder with kilosort results (Type Z: to change to engram drive. then cd path. )
%   phy template-gui params.py


%% Get waveforms for all sessions in table T
only_good = 0; % Process all for now, select good later
for ii = 1:height(T)
    
    % Find the binary directory
    root_dir = fullfile(data_dir_local, T.filename{ii});
    raw_dir_temp = dir(fullfile(root_dir, 'raw*'));
    raw_dir = fullfile(raw_dir_temp.folder, raw_dir_temp.name);
    ks_dir = ksDirFun(root_dir);
    
    % Check if already done with the latest sorting unit labels
    if exist(fullfile(ks_dir, 'waveformStruct.mat'),'file')
        wvStruct = getfield(load(fullfile(ks_dir, 'waveformStruct.mat')),'wvStruct');
        
        % Get the latest phy labels
        [unit_ID,cluster_labels] = getPhyClusterLabels(ks_dir);
        
        % Check if everything is already identical
        if length(unit_ID)==length(wvStruct.goodIDs) && all(unit_ID(:)==wvStruct.goodIDs(:)) && all(strcmp(wvStruct.goodLabels,cluster_labels))
            continue;
            
            % Check if units are identical, but labels (good/mua/noise) need to be updated
        elseif length(unit_ID)==length(wvStruct.goodIDs) && all(unit_ID(:)==wvStruct.goodIDs(:)) && ~all(strcmp(wvStruct.goodLabels,cluster_labels))
            % Re-save in case any updates to KS labels
            wvStruct.goodLabels = cluster_labels;
            save(fullfile(ks_dir,'waveformStruct.mat'),'wvStruct')
            continue;
        end
        
    end
    
    % Save waveforms for cell type clustering
    disp(T(ii,:))
    wvStruct = getSessionWaveforms(raw_dir, ks_dir, only_good);
    save(fullfile(ks_dir, 'waveformStruct.mat'), 'wvStruct','-v7.3')
    save(fullfile(ksDirFun(fullfile(data_dir_remote, T.filename{ii})), 'waveformStruct.mat'), 'wvStruct','-v7.3')
    
end

%% Generate GMM based on all curated sessions - rerun after sorting new sessions
gmm_name = 'latestGMM.mat';
if ~exist(fullfile(results_dir, gmm_name),'file') || overwrite_gmm
    T_load = T(strcmp(T.manually_sorted,'yes'),:);
    option_only_good = 1;
    n_clusters = 2;
    gm = GMM_make(cellfun(@(x) fullfile(data_dir_local, x), T_load.filename,'Uni',0), option_only_good, n_clusters);
    save(fullfile(results_dir, gmm_name),'gm','T')
end


%% Apply GMM results to all sessions
gm = getfield(load(fullfile(results_dir, gmm_name)),'gm');
option_only_good = 0;
plot_on = 0;
for ii = 1:height(T)
    ks_dir = ksDirFun(fullfile(data_dir_local, T.filename{ii}));
    if ~exist(fullfile(ks_dir,'gmm_result.mat'),'file') || overwrite_gmm
        disp(T.filename{ii})
        [idx,labels] = GMM_apply(gm, fullfile(data_dir_local, T.filename{ii}), option_only_good, plot_on);
        save(fullfile(ks_dir,'gmm_result.mat'),'idx','labels')
        
        wvStruct = getfield(load(fullfile(ks_dir, 'waveformStruct.mat')), 'wvStruct');
        wvStruct.typeLabels = labels;
        save(fullfile(ks_dir, 'waveformStruct.mat'),'wvStruct')
    end
end

%% Apply basic sorting quality metrics to uncurated sessions
% Make sure to record if manually sorted in the excel table!!!
for ii = 1:height(T)
    ks_dir = ksDirFun(fullfile(data_dir_local, T.filename{ii}));
    if strcmp(T.manually_sorted{ii},'yes'); disp('Skipping, manually sorted');
        continue;
    end
    applyQualityMetrics(ks_dir);
end

%% Backup kilosort results
exceptions = {'params.py','phy.log','.phy'};
runmode = 0;
for ii = 1:height(T)
    fprintf('\n%s\n',T.filename{ii})
    s1 = ksDirFun(fullfile(data_dir_local, T.filename{ii}));
    s2 = ksDirFun(fullfile(data_dir_remote, T.filename{ii}));
    dirbackup(s1, s2, runmode, exceptions)
end

%% Optional plotting of # of cells and E/I ratio
plotBasicResults(T, data_dir_local, ksDirFun)
