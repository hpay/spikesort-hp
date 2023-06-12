%% Generate list of directories to run kilosort on
%
% Instructions: Update RECORDING_DEPTH_CHICK.xlsx then run this script
raw_dir_temp = which(mfilename);
dropbox_folder = raw_dir_temp(1:strfind(raw_dir_temp,'Dropbox')+6);
% T = readtable(fullfile(dropbox_folder,'alab\Analysis\RECORDING_DEPTH_CHICK.xlsx'));
T = readtable(fullfile(dropbox_folder,'alab\Code\project2\data\RECORDING_DEPTH_CHICK.xlsx'));
T = T(T.include  & (strcmp(T.task,'X_gaze')|strcmp(T.task,'X_no_gaze')|strcmp(T.task,'X_fix')) ,:);
fs = 3e4;

% Overwrite kilosort output?
overwrite = 0;

% Data folder
data_dir_remote = 'Z:\Hannah\Ephys\Project2';
data_dir_local = 'D:\data'; 

% Save dir
ks_dir_fun = @(root_dir) fullfile(root_dir,'kilosort2_output');

% Results folder (GMM etc)
addpath('..\results')

%% set up code folders
% Copy spike sorting results and raw data for post-processing and analysis. NOT backed up - make sure to copy anything saved here back to server or have code to reproduce!
scripts_dir = cd;
spikesort_hp_dir = fileparts(scripts_dir);
code_dir = fileparts(spikesort_hp_dir);
addpath(genpath(fullfile(spikesort_hp_dir,'src')))
warning("off","parallel:gpu:device:DeviceDeprecated");
addpath(genpath(fullfile(code_dir, 'kilosort-2.0'))) % Add kilosort2 directory

%% Backup from local SSD to server
exceptions = {'kilosort2_output'}; % Back this up specifically later
runmode = 0;
for ii = 1:height(T)
    fprintf('\n%s\n',T.filename{ii})
    root_dir_local = fullfile(data_dir_local, T.filename{ii});  % e.g. D:\data\HC11_230129 
    root_dir_remote = fullfile(data_dir_remote, T.filename{ii}); % e.g. Z:\Hannah\ephys\HC11_230129 
    mkdir(root_dir_remote)
    dirbackup(root_dir_local, root_dir_remote, runmode, exceptions)
end

%% Run Kilosort2 
for ii = 1:height(T)
    % Root dir on server (Z:\Hannah\ephys\HC11_230129 etc)
    root_dir_remote = fullfile(data_dir_remote, T.filename{ii});
    root_dir_local = fullfile(data_dir_local, T.filename{ii});

    % Find the binary directory
    raw_dir_temp = dir(fullfile(root_dir_local, 'raw*'));
    raw_dir = fullfile(raw_dir_temp.folder, raw_dir_temp.name);
    
    % Create save directory
    ks_dir = ks_dir_fun(root_dir_local);
    
    % If it exists already, either overwrite or skip
    if ~isempty(dir(fullfile(ks_dir,'*.npy')))
        if overwrite
            fprintf('Overwriting %s\n', T.filename{ii})
            try
                rmdir(ks_dir,'s')
            catch
            end
        else
            fprintf('Skipping %s\n', T.filename{ii})
            continue;
        end
    end
    mkdir(ks_dir)
    disp(T(ii,:))
    
    % set up options
    ops = [];
    ops.badChannels = eval(T.bad_chan{ii});
    ops.chanMap = [T.probe_chanmap{ii} '.mat']; % Make sure config folder is on the search path
    ops = hp_config(ops);
    
    % Run and save results. Note: ops saved in rez.mat
    fprintf('\n Running Kilosort on directory %s \n', raw_dir),
    run_single_kilosort(raw_dir, ks_dir, ops);       
    
    % Copy kilosort output to the server
    save_dir_remote = ks_dir_fun(root_dir_remote);
    copyfile(ks_dir, save_dir_remote)   
    
    % Edit params.py to reflect new amplifier.dat location
    fid = fopen(fullfile(ks_dir,'params.py'), 'r'); % Read params.py from local copy
    a = {};
    while ~feof(fid); a{end+1} = fgets(fid);  pause(.001);  end
    fclose(fid);
    a{1} = sprintf("dat_path = '%s'\n", strrep(fullfile(root_dir_remote, raw_dir_temp.name, 'amplifier.dat'),'\','/'));
    fid = fopen(fullfile(save_dir_remote,'params.py'), 'w'); % Write correct params.py to remote server
    for jj = 1:length(a); fprintf(fid,a{jj});      end
    fclose(fid);
    fprintf('Results copied to server!\n')
    
end

%% Manually label results in Phy! 
% disp('Manually sort now!')
% keyboard


%% Get waveforms for all sessions in table T
only_good = 0; % Process all for now, select good later
for ii = 1:height(T)
    
    % Find the binary directory
    root_dir = fullfile(data_dir_local, T.filename{ii});
    raw_dir_temp = dir(fullfile(root_dir, 'raw*'));
    raw_dir = fullfile(raw_dir_temp.folder, raw_dir_temp.name);
    ks_dir = ks_dir_fun(root_dir);
    
    % Check if already done with the latest sorting unit labels
    if exist(fullfile(ks_dir, 'waveformStruct.mat'),'file')
        wvStruct = getfield(load(fullfile(ks_dir, 'waveformStruct.mat')),'wvStruct');

        % Get the latest phy labels
        [unit_ID,cluster_labels] = get_phy_cluster_labels(ks_dir);
                
        % Check if everything is already identical
        if length(unit_ID)==length(wvStruct.goodIDs) && all(unit_ID(:)==wvStruct.goodIDs(:)) && all(strcmp(wvStruct.goodLabels,cluster_labels))
            continue;
            
        % Check if units are identical, but labels (good/mua/noise) need to be updated
        elseif length(unit_ID)==length(wvStruct.goodIDs) && all(unit_ID(:)==wvStruct.goodIDs(:)) && ~all(strcmp(wvStruct.goodLabels,cluster_labels))
            % Resave in case any updates to KS labels
            wvStruct.goodLabels = cluster_labels;
            save(fullfile(ks_dir,'waveformStruct.mat'),'wvStruct')
            continue;
        end        
        
    end
    
    % Save waveforms for cell type clustering
    disp(T(ii,:))
    wvStruct = getSessionWaveforms(raw_dir, ks_dir, only_good);
    save(fullfile(ks_dir, 'waveformStruct.mat'), 'wvStruct','-v7.3')
    save(fullfile(ks_dir_fun(fullfile(data_dir_remote, T.filename{ii})), 'waveformStruct.mat'), 'wvStruct','-v7.3')
    
end

%% Generate GMM based on all curated sessions - rerun after sorting new sessions
%{
T_load = T(strcmp(T.manually_sorted,'yes'),:);
option_only_good = 1;
n_clusters = 2;
gm = GMM_make(cellfun(@(x) fullfile(data_dir_local, x), T_load.filename,'Uni',0), option_only_good, n_clusters);
save('..\results\latestGMM_cellType.mat','gm','T')
%}

%% Apply GMM results to all sessions
load('..\results\latestGMM_cellType.mat','gm')
option_only_good = 0;
plot_on = 0;
for ii = 1:height(T)
    disp(T.filename{ii})
    ks_dir = ks_dir_fun(fullfile(data_dir_local, T.filename{ii}));
    [idx,labels] = GMM_apply(gm, fullfile(data_dir_local, T.filename{ii}), option_only_good, plot_on);
    save(fullfile(ks_dir,'gmm_result.mat'),'idx','labels')
end

%% Apply basic sorting quality metrics
max_contam_good = 0.10;
max_contam_mua = 1;
min_spikes_good = 50; % Label units with less than 50 spikes as mua
min_spikes_mua = 20; % Label units with less than 20 spikes as noise 

for ii = 1:height(T)
    % Find the KS directory
    ks_dir = ks_dir_fun(fullfile(data_dir_local, T.filename{ii}));
        
    % Check if manually curated - make sure to record in the excel table!!!
    disp(T.filename{ii})
    
    % Get the contam percents calculated in getSessionWaveforms
    wvStruct = getfield(load(fullfile(ks_dir, 'waveformStruct.mat')), 'wvStruct');
    gm_result = load(fullfile(ks_dir,'gmm_result.mat'));
    wvStruct.typeLabels = gm_result.labels;
    save(fullfile(ks_dir, 'waveformStruct.mat'),'wvStruct')
    
    if strcmp(T.manually_sorted{ii},'yes'); disp('Skipping, manually sorted');
        continue;
    end
    
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
    
    
    % Label units according to GMM
    mask_good = mask_good & ismember(gm_result.labels,{'E','I'});
    mask_noise = mask_noise | ismember(gm_result.labels,{'extreme_outlier','unknown'});
    mask_mua = ~mask_noise & ~mask_good;

    % TODO: Label units according to firing rate stability
    
    if nnz(mask_mua)+nnz(mask_good)+nnz(mask_noise) ~= length(mask_noise); error('check masks'); end
    
%     figure('Pos',[209 266 1065 300]); subplot(1,3,1); plot(wvStruct.mxWF(mask_good,:)','g'); subplot(1,3,2);  plot(wvStruct.mxWF(mask_mua,:)','k');  subplot(1,3,3); plot(wvStruct.mxWF(mask_noise,:)','r');
%     linkaxes; ylim([-1.2 .6]*1e3); 
%     title(T.filename{ii},'Interp','none');
    
    % Save new "good", "mua", "noise" labels
    KSLabel = cell(length(cluster_id),1);
    KSLabel(mask_good) = deal({'good'});
    KSLabel(mask_mua) = deal({'mua'});
    KSLabel(mask_noise) = deal({'noise'});
    fprintf('%s: ngood new %i, ngood orig %i\n', T.filename{ii}, nnz(mask_good), nnz(strcmp(T_cluster_KSLabel_old.KSLabel,'good')))

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
end

%% Plot fraction E v. depth - good+ mua
nE = []; nI = [];
 for ii = 1:height(T)
     ks_dir = ks_dir_fun(fullfile(data_dir_local, T.filename{ii}));
     wvStruct = getfield(load(fullfile(ks_dir, 'waveformStruct.mat')), 'wvStruct');     
%      nE(ii) = nnz(strcmp(wvStruct.goodLabels,'good') & strcmp(wvStruct.typeLabels,'E'));
%      nI(ii) = nnz(strcmp(wvStruct.goodLabels,'good') & strcmp(wvStruct.typeLabels,'I'));
     nE(ii) = nnz(ismember(wvStruct.goodLabels,{'good','mua'}) & strcmp(wvStruct.typeLabels,'E'));
     nI(ii) = nnz(ismember(wvStruct.goodLabels,{'good','mua'}) & strcmp(wvStruct.typeLabels,'I'));
 end
 T.nE = nE(:); T.nI= nI(:);
T.ratioE = T.nE./(T.nE+T.nI);
disp(T(:,{'filename','depth','nE','nI','ratioE'}))

figure;
cs = turbo(7);
birds = unique(T.bird);
hs = gobjects(length(birds),1);
for ii = 1:length(birds)
    mask = strcmp(T.bird,birds{ii});
    hs(ii) = plot(T.depth(mask), T.ratioE(mask),'o','Color',...
        cs(strcmp(birds, birds{ii}),:),'LineWidth',2); hold on
%         hs(ii) = plot(T.depth(mask), T.nE(mask)+T.nI(mask),'o','Color',c(strcmp(birds, birds{ii}),:),'LineWidth',2); hold on
end
xlabel('Depth (tip, mm)')
ylabel('Fraction of cells E')
% ylabel('Total cells E+I')
grid
legend(birds)
ylim([0 1])

%% Plot the number of "good" cells (manual and KS labelled) v. days since implant
disp(T(:,{'filename','depth','manually_sorted','probe_chanmap','bad_chan','notes'}))
disp('Counting total cells')
ngood_final= [];
ngood_orig = [];
for ii = 1:height(T)
    % Find the binary directory
    root_dir = fullfile(data_dir_local, T.filename{ii});
    raw_dir_temp = dir(fullfile(root_dir, 'raw*'));
    raw_dir = fullfile(raw_dir_temp.folder, raw_dir_temp.name);
    ks_dir = ks_dir_fun(root_dir);

    % Latest labels saved by Phy
    [unit_ID,cluster_labels] = get_phy_cluster_labels(ks_dir);
    ngood_final(ii) = nnz(strcmp(cluster_labels,'good'));

    % Original good labels (stored in 'rez.mat')
    good_orig = getfield(load(fullfile(ks_dir,'rez.mat'),'good'),'good');
    ngood_orig(ii) = nnz(good_orig);
    
end

% Plot the number of "good" cells per session
days_since_surgery = days(T.date-T.surgery);
birds = unique(T.bird);
figure;
hold on
% cs = hsv(length(birds)+1);
hs = gobjects(length(birds),2);
for jj = 1:length(birds)
    mask =strcmp(T.bird, birds{jj});
    mask_manual = strcmp(T.manually_sorted,'yes');
    hs(jj,1) = plot(days_since_surgery(mask), ngood_orig(mask),'-o','Color',cs(jj,:),'MarkerFaceColor','w','MarkerEdgeColor',cs(jj,:));
    if any(mask_manual&mask)
    plot(days_since_surgery(mask&mask_manual), ngood_final(mask&mask_manual),'o','Color',cs(jj,:),'MarkerFaceColor',cs(jj,:),'MarkerEdgeColor',cs(jj,:));
    end
    hold on;
end
xlabel('Days since surgery'); ylabel('# "good" cells')
ylim([0 max(ngood_orig)+2])
grid on
legend(hs(:,1),birds)

%% Backup kilosort results
exceptions = {'params.py','phy.log','.phy'};
runmode = 0;
for ii = 1:height(T)
    fprintf('\n%s\n',T.filename{ii})
    s1 = ks_dir_fun(fullfile(data_dir_local, T.filename{ii}));
    s2 = ks_dir_fun(fullfile(data_dir_remote, T.filename{ii}));
    dirbackup(s1, s2, runmode, exceptions)
end