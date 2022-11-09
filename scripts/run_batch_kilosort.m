%% Generate list of directories to run kilosort on
T = readtable('D:\hannah\Dropbox\alab\Analysis\RECORDING_DEPTH_CHICK.xlsx');
T = T(~T.exclude,:);
fs = 3e4;

% Overwrite kilosort output?
overwrite = 0;

% Data folder
data_dir = 'Z:\Hannah\Ephys\Project2';

% Save dir
save_dir_fun = @(root_dir) fullfile(root_dir,'kilosort2_output');

%% set up code folders
scripts_dir = cd;
spikesort_hp_dir = fileparts(scripts_dir);
code_dir = fileparts(spikesort_hp_dir);
addpath(genpath(fullfile(spikesort_hp_dir,'src'))) 
warning("off","parallel:gpu:device:DeviceDeprecated");
addpath(genpath(fullfile(code_dir, 'kilosort-2.0'))) % Add kilosort2 directory

%% Run all sessions in table T
for ii = 1:height(T)
    disp(T(ii,:))
    
    % Find the binary directory
    root_dir = fullfile(data_dir, T.filename{ii});
    temp = dir(fullfile(data_dir, T.filename{ii}, 'raw*'));
    raw_dir = fullfile(temp.folder, temp.name);
    
    % Create save directory
    save_dir = save_dir_fun(root_dir);
    mkdir(save_dir);
        
    % If it exists already, either overwrite or skip
    if exist(save_dir,'dir')
        if overwrite
            disp('overwriting')
            try
                rmdir(save_dir,'s')
            catch
            end
        else
            disp('skipping')
            continue;
        end
    end
    mkdir(save_dir)
    
    % set up options
    ops = [];
    ops.badChannels = eval(T.bad_chan{ii});
    ops.chanMap = [T.probe_chanmap{ii} '.mat']; % Make sure config folder is on the search path
    ops = hp_config(ops);
    
    % Run and save results. Note: ops saved in rez.mat
    fprintf('\n Running Kilosort on directory %s \n', raw_dir),
    run_single_kilosort(raw_dir, save_dir, ops);    
    
end

%% Manually label results in Phy! Then run waveforms


%% Get waveforms for all sessions in table T
T = T(strcmp(T.manually_sorted,'yes'),:)

for ii = 1%:height(T)
    
        % Find the binary directory
    disp(T(ii,:))
    root_dir = fullfile(data_dir, T.filename{ii});
    temp = dir(fullfile(data_dir, T.filename{ii}, 'raw*'));
    raw_dir = fullfile(temp.folder, temp.name);
        save_dir = save_dir_fun(root_dir);

    % Save waveforms for celltype clustering
    wvStruct = getSessionWaveforms(raw_dir, save_dir, 0);
    save(fullfile(save_dir, 'waveformStruct.mat'), 'wvStruct')
   
%     wvStruct = getfield(load(fullfile(save_dir,'waveformStruct.mat')),'wvStruct');
% [unit_ID, cluster_labels] = get_phy_cluster_labels(save_dir); % 0 = noise, 1 = 'mua', 2 = 'good'

end


%% Load all spikes sorted and plot waveforms and clusters
T = T(strcmp(T.manually_sorted,'yes'),:)
ks_dir_fun = @(root_dir) fullfile(root_dir,'kilosort2_output');
option_only_good = 1; % Only include cells I've marked as good

clear S
for ii = 1:height(T)
    save_dir = ks_dir_fun(fullfile(data_dir, T.filename{ii}));      
    S{ii} =  importKilosort(save_dir, fs, option_only_good);        
end

%% Plot the number of "good" cells per session
ncells = cellfun(@length,S);
days_since_surgery = days(T.date-T.surgery);
birds = unique(T.bird);
figure;
hold on
for jj = 1:length(birds)
    mask =strcmp(T.bird, birds{jj});
    plot(days_since_surgery(mask), ncells(mask),'-o')
    hold on;
end
xlabel('Days since surgery'); ylabel('# "good" cells in Kilosort2 (SC params)')
ylim([0 max(ncells)+2])
xlim([0 16])
grid on 
legend(birds)

%% Plot waveforms
figure
for ii = 1:height(T)
    waves = cell2mat({S{ii}.waveform}')';
    for jj = 1:size(waves,2)
        [maxval, maxind] = max(abs(waves(:,jj)));
        plot(scaledata(waves(:,jj),0, maxval,0, 1),'k'); hold on
    end
end

%% Generate GMM for all sorted cells
gm = makeGMM_celltypes(cellfun(@(x) fullfile(data_dir, x), T.filename,'Uni',0))
save('..\results\latestGMM_cellType.mat','gm','T')



