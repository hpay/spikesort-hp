%% Generate list of directories to run kilosort on

raw_dir_temp = which(mfilename);
dropbox_folder = raw_dir_temp(1:strfind(raw_dir_temp,'Dropbox')+6);
T = readtable(fullfile(dropbox_folder,'alab\Analysis\RECORDING_DEPTH_CHICK.xlsx'));
T = T(T.include  & (strcmp(T.task,'X_gaze')|strcmp(T.task,'X_no_gaze')) ,:);
fs = 3e4;

% Overwrite kilosort output?
overwrite = 0;

% Data folder
data_dir_remote = 'Z:\Hannah\Ephys\Project2';
data_dir_local = 'D:\temp'; 

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

%% Run all sessions in table T
for ii = 1:height(T)
    
    % Root dir on server (Z:\Hannah\ephys\HC11_230129 etc)
    root_dir_remote = fullfile(data_dir_remote, T.filename{ii});
    root_dir_local = fullfile(data_dir_local, T.filename{ii});

    % Copy everything to local data SSD
    if ~exist(root_dir_local,'dir') && exist(root_dir_remote,'dir')
        copyfile(root_dir_remote, root_dir_local)
    end
        
end


for ii = 1:height(T)
    % Root dir on server (Z:\Hannah\ephys\HC11_230129 etc)
    root_dir_remote = fullfile(data_dir_remote, T.filename{ii});
    root_dir_local = fullfile(data_dir_local, T.filename{ii});

    % Find the binary directory
    raw_dir_temp = dir(fullfile(root_dir_local, 'raw*'));
    raw_dir = fullfile(raw_dir_temp.folder, raw_dir_temp.name);
    
    % Create save directory
    save_dir = ks_dir_fun(root_dir_local);
    
    % If it exists already, either overwrite or skip
    if ~isempty(dir(fullfile(save_dir,'*.npy')))
        if overwrite
            fprintf('Overwriting %s\n', T.filename{ii})
            try
                rmdir(save_dir,'s')
            catch
            end
        else
            fprintf('Skipping %s\n', T.filename{ii})
            continue;
        end
    end
    mkdir(save_dir)
    disp(T(ii,:))
    
    % set up options
    ops = [];
    ops.badChannels = eval(T.bad_chan{ii});
    ops.chanMap = [T.probe_chanmap{ii} '.mat']; % Make sure config folder is on the search path
    ops = hp_config(ops);
    
    % Run and save results. Note: ops saved in rez.mat
    fprintf('\n Running Kilosort on directory %s \n', raw_dir),
    run_single_kilosort(raw_dir, save_dir, ops);
        
    
    % Apply my custom prescreening and apply labels to 
    
    % Copy kilosort output to the server
    save_dir_remote = ks_dir_fun(root_dir_remote);
    copyfile(save_dir, save_dir_remote)   
    
    % Edit params.py to reflect new amplifier.dat location
    fid = fopen(fullfile(save_dir,'params.py'), 'r'); % Read params.py from local copy
    a = {};
    while ~feof(fid); a{end+1} = fgets(fid);  pause(.001);  end
    fclose(fid);
    a{1} = sprintf("dat_path = '%s'\n", strrep(fullfile(root_dir_remote, raw_dir_temp.name, 'amplifier.dat'),'\','/'));
    fid = fopen(fullfile(save_dir_remote,'params.py'), 'w'); % Write correct params.py to remote server
    for jj = 1:length(a); fprintf(fid,a{jj});      end
    fclose(fid);
    
    fprintf('Results copied to server!\n')
    
end

%% Manually label results in Phy! Then run waveforms
% disp('Manually sort now!')
% keyboard
data_dir_process = data_dir_local

% T = T([1:20],:)
%% Plot the number of "good" cells (manual and KS labelled)
disp(T(:,{'filename','depth','manually_sorted','probe_chanmap','bad_chan','notes'}))
disp('Counting total cells')
ngood_final= [];
ngood_orig = [];
for ii = 1:height(T)
    % Find the binary directory
    root_dir = fullfile(data_dir_process, T.filename{ii});
    raw_dir_temp = dir(fullfile(root_dir, 'raw*'));
    raw_dir = fullfile(raw_dir_temp.folder, raw_dir_temp.name);
    save_dir = ks_dir_fun(root_dir);

    % Latest labels saved by Phy
    [cIDs,cluster_labels] = get_phy_cluster_labels(save_dir);
    ngood_final(ii) = nnz(strcmp(cluster_labels,'good'));

    % Original good labels (stored in 'rez.mat')
    good_orig = getfield(load(fullfile(save_dir,'rez.mat'),'good'),'good');
    ngood_orig(ii) = nnz(good_orig);
    
end

% Plot the number of "good" cells per session
days_since_surgery = days(T.date-T.surgery);
birds = unique(T.bird);
figure;
hold on
cs = hsv(length(birds)+1);
hs = gobjects(length(birds),2);
for jj = 1:length(birds)
    mask =strcmp(T.bird, birds{jj});
    hs(jj,1) = plot(days_since_surgery(mask), ngood_final(mask),'o','Color',cs(jj,:),'MarkerFaceColor',cs(jj,:),'MarkerEdgeColor',cs(jj,:));
    hs(jj,2) = plot(days_since_surgery(mask), ngood_orig(mask),'-o','Color',cs(jj,:),'MarkerFaceColor','w','MarkerEdgeColor',cs(jj,:));
    
    hold on;
end
xlabel('Days since surgery'); ylabel('# "good" cells')
ylim([0 max(ngood_orig)+2])
% xlim([0 16])
grid on
legend(hs(:,1),birds)


%% Get waveforms for all sessions in table T
for ii = 1:height(T)
    
    % Find the binary directory
    root_dir = fullfile(data_dir_process, T.filename{ii});
    raw_dir_temp = dir(fullfile(root_dir, 'raw*'));
    raw_dir = fullfile(raw_dir_temp.folder, raw_dir_temp.name);
    save_dir = ks_dir_fun(root_dir);
    
    % Check if already done with the latest sorting unit labels
    if exist(fullfile(save_dir, 'waveformStruct.mat'),'file')
        wvStruct =getfield(load(fullfile(save_dir, 'waveformStruct.mat')),'wvStruct');
        [cIDs,cluster_labels] = get_phy_cluster_labels(save_dir);
        if length(wvStruct)==1 && length(wvStruct.goodLabels)>1 && length(cIDs)==length(wvStruct.goodIDs) && all(cIDs==wvStruct.goodIDs);
            continue;
        end
        
    end
    % Save waveforms for celltype clustering
    disp(T(ii,:))
    only_good = 0; % Process all for now, select good later
    wvStruct = getSessionWaveforms(raw_dir, save_dir, only_good);
    save(fullfile(save_dir, 'waveformStruct.mat'), 'wvStruct','-v7.3')
    save(fullfile(ks_dir_fun(fullfile(data_dir_remote, T.filename{ii})), 'waveformStruct.mat'), 'wvStruct','-v7.3')
    
end


%% Load all spikes sorted and plot waveforms and clusters
T_load = T(strcmp(T.manually_sorted,'yes'),:)

option_only_good = 1; % Only include cells marked as good

% clear S
for ii = 1:height(T_load)
    save_dir = ks_dir_fun(fullfile(data_dir_local, T_load.filename{ii}));
    S{ii} =  importKilosort(save_dir, fs, option_only_good);
end


% Plot waveforms
figure
c = jet(height(T_load));
hs = gobjects(height(T_load),1);
for ii = 1:height(T_load)
    waves = cell2mat({S{ii}.waveform}')';
    for jj = 1:size(waves,2)
        [maxval, maxind] = max(abs(waves(:,jj)));
        hs(ii) = plot(scaledata(waves(:,jj),0, maxval,0, 1),'Color',c(ii,:)); hold on
    end
end
legend(hs,arrayfun(@(x) num2str(x), 1:height(T_load),'Uni',0))

%% Generate GMM for all sorted cells
Tmanual = T(strcmp(T.manually_sorted,'yes'),:)

option_only_good = 1;
gm = makeGMM_celltypes(cellfun(@(x) fullfile(data_dir_local, x), Tmanual.filename,'Uni',0), option_only_good);
save('..\results\latestGMM_cellType.mat','gm','T')

%% Save cell identity results
option_only_good = 0;
for ii = 1:height(T)
    save_dir_local = ks_dir_fun(fullfile(data_dir_local, T.filename{ii}));
    save_dir = ks_dir_fun(fullfile(data_dir_remote, T.filename{ii}));
    [idx,nlogL,P,logpdf, goodIDs] = classifyUnitsGMM(fullfile(data_dir_local, T.filename{ii}), option_only_good);
    save(fullfile(save_dir,'gmm_result.mat'),'idx','nlogL','P','logpdf', 'goodIDs')
    save(fullfile(save_dir,'gmm_result.mat'),'idx','nlogL','P','logpdf', 'goodIDs')
end