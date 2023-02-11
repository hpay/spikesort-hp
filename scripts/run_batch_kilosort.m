%% Generate list of directories to run kilosort on

temp = which(mfilename);
dropbox_folder = temp(1:strfind(temp,'Dropbox')+6);
T = readtable(fullfile(dropbox_folder,'alab\Analysis\RECORDING_DEPTH_CHICK.xlsx'));
T = T(~T.exclude,:);
fs = 3e4;

% Overwrite kilosort output?
overwrite = 0;

% Data folder
data_dir = 'Z:\Hannah\Ephys\Project2';

% Save dir
save_dir_fun = @(root_dir) fullfile(root_dir,'kilosort2_output');

% Results folder (GMM etc)
addpath('..\results')

%% set up code folders
scripts_dir = cd;
spikesort_hp_dir = fileparts(scripts_dir);
code_dir = fileparts(spikesort_hp_dir);
addpath(genpath(fullfile(spikesort_hp_dir,'src')))
warning("off","parallel:gpu:device:DeviceDeprecated");
addpath(genpath(fullfile(code_dir, 'kilosort-2.0'))) % Add kilosort2 directory

%% Run all sessions in table T
for ii = 1:height(T)
    
    % Find the binary directory
    root_dir = fullfile(data_dir, T.filename{ii});
    temp = dir(fullfile(data_dir, T.filename{ii}, 'raw*'));
    raw_dir = fullfile(temp.folder, temp.name);
    
    % Create save directory
    save_dir = save_dir_fun(root_dir);
    
    % If it exists already, either overwrite or skip
    if exist(save_dir,'dir') && ~isempty(dir(fullfile(save_dir,'*.npy')))
        if overwrite
            
            fprint('Overwriting %s\n', T.filename{ii})
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
    
    
    % OPTIONAL: copy raw data and kilosort output to a temp folder on a
    % local drive
    temp_folder = 'D:\temp';
    mkdir(fullfile(temp_folder, T.filename{ii}))
    copyfile(raw_dir, fullfile(temp_folder, T.filename{ii}))
    copyfile(save_dir_fun(fullfile(data_dir, T.filename{ii})), ...
        save_dir_fun(fullfile(temp_folder, T.filename{ii})))
    fprintf('Results copied to %s - delete raw files after manually inspecting!\n', temp_folder)
    
end

%% Manually label results in Phy! Then run waveforms
% disp('Manually sort now!')
% keyboard



%% Plot the number of "good" cells (manual and KS labelled(
disp(T(:,{'filename','depth','manually_sorted','probe_chanmap','bad_chan','notes'}))
disp('Counting total cells')
ngood_final= [];
ngood_orig = [];
for ii = 1:height(T)
    % Find the binary directory
    root_dir = fullfile(data_dir, T.filename{ii});
    temp = dir(fullfile(data_dir, T.filename{ii}, 'raw*'));
    raw_dir = fullfile(temp.folder, temp.name);
    save_dir = save_dir_fun(root_dir);

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
% **** TODO: fix this! why is wvStruct repeated nUnit times with the same
% info in each!
for ii = 1:height(T)
    
    % Find the binary directory
    root_dir = fullfile(data_dir, T.filename{ii});
    temp = dir(fullfile(data_dir, T.filename{ii}, 'raw*'));
    raw_dir = fullfile(temp.folder, temp.name);
    save_dir = save_dir_fun(root_dir);
    
    % Check if already done with the latest sorting unit labels
    if exist(fullfile(save_dir, 'waveformStruct.mat'),'file')
        wvStruct =getfield(load(fullfile(save_dir, 'waveformStruct.mat')),'wvStruct');
        [cIDs,cluster_labels] = get_phy_cluster_labels(save_dir);
        if length(cIDs)==length(wvStruct.goodIDs) && all(cIDs==wvStruct.goodIDs);
            continue;
        end
        
    end
    % Save waveforms for celltype clustering
    disp(T(ii,:))
    only_good = 0; % Process all for now, select good later
    wvStruct = getSessionWaveforms(raw_dir, save_dir, only_good);
    save(fullfile(save_dir, 'waveformStruct.mat'), 'wvStruct','-v7.3')
    
end





%% Load all spikes sorted and plot waveforms and clusters
T = T(strcmp(T.manually_sorted,'yes'),:)

ks_dir_fun = @(root_dir) fullfile(root_dir,'kilosort2_output');
option_only_good = 1; % Only include cells marked as good

% clear S
for ii = 1:height(T)
    save_dir = ks_dir_fun(fullfile(data_dir, T.filename{ii}));
    S{ii} =  importKilosort(save_dir, fs, option_only_good);
end


%% Plot waveforms
figure
c = jet(height(T));
hs = gobjects(height(T),1);
for ii = 1:height(T)
    waves = cell2mat({S{ii}.waveform}')';
    for jj = 1:size(waves,2)
        [maxval, maxind] = max(abs(waves(:,jj)));
        hs(ii) = plot(scaledata(waves(:,jj),0, maxval,0, 1),'Color',c(ii,:)); hold on
    end
end
legend(hs,arrayfun(@(x) num2str(x), 1:height(T),'Uni',0))

%% Generate GMM for all sorted cells
T = T(strcmp(T.manually_sorted,'yes'),:)

option_only_good = 1;
gm = makeGMM_celltypes(cellfun(@(x) fullfile(data_dir, x), T.filename,'Uni',0), option_only_good);
save('..\results\latestGMM_cellType.mat','gm','T')

%% Save cell identity results
option_only_good = 0;
for ii = 1:height(T)
    save_dir = ks_dir_fun(fullfile(data_dir, T.filename{ii}));
    [idx,nlogL,P,logpdf, goodIDs] = classifyUnitsGMM(fullfile(data_dir, T.filename{ii}), option_only_good);
    save(fullfile(save_dir,'gmm_result.mat'),'idx','nlogL','P','logpdf', 'goodIDs')
end