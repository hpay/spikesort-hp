%% Generate list of directories to run kilosort on
% birdID = 'HC05'; % AMB111
T = readtable('D:\hannah\Dropbox\alab\Analysis\RECORDING_DEPTH_CHICK.xlsx');
T = T(~T.exclude,:);
% T = T(strcmp(T.bird, birdID),:)
fs = 3e4;

% Overwrite kilosort output?
overwrite = 1;

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

%% Add kilosort2 directory
addpath(genpath(fullfile(code_dir, 'kilosort-2.0')))

%% Run a single session 
%{
% raw_dir = 'Z:\Hannah\ephys\project2\HC05_220819\raw_220819_125309'
raw_dir = 'Z:\Hannah\ephys\project2\HC09_221104\raw_221104_121402';
ops = [];
ops.chanMap = 'H6.mat';
% ops.chanMap = 'H6_shankA.mat';
ops = hp_config(ops);

save_dir = save_dir_fun(fileparts(raw_dir)); 
save_dir= [save_dir '_minfr_goodchannels1_50']
mkdir(save_dir)
fprintf('\n Running Kilosort on directory %s \n', raw_dir)
t = tic;
run_single_kilosort(raw_dir,save_dir,ops)
toc(t)

% Save waveforms for celltype clustering
wvStruct = getSessionWaveforms(raw_dir, save_dir, 0);
save(fullfile(save_dir, 'waveformStruct.mat'), 'wvStruct')

% Load spikes into my dat structure in matlab
option_only_good= 1;
S = importKilo2(save_dir, ops.fs, option_only_good);
disp(S)

% Run phy to inspect it! (TODO)
% In Anaconda prompt:
%   conda activate phy_env
% Navigate to folder with kilosort results (Type Z: to change to engram drive. then cd path. )
%   phy template-gui params.py
%}

%% Run all sessions in table T
for ii = height(T):-1:1
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
            catch;
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
    
    % Save waveforms for celltype clustering
    wvStruct = getSessionWaveforms(raw_dir, save_dir, 0);
    save(fullfile(save_dir, 'waveformStruct.mat'), 'wvStruct')
    %     catch msg
    %         disp(msg)
    %     end
end


%% Collect # of good cells in each?
clear S
for ii = 1:height(T)

    % Find the save directory
    save_dir = save_dir_fun(fullfile(data_dir, T.filename{ii}));      
    option_only_good = 1;
    S{ii} =  importKilo2(save_dir, fs, option_only_good);
        
end
ncells = cellfun(@length,S);
days_since_surgery = days(T.date-datetime(2022,08,09));

figure;
hold on
plot(days_since_surgery, ncells,'-ob')
xlabel('Days since surgery'); ylabel('# "good" cells in Kilosort2 (SC params)')
ylim([0 max(ncells)+2])
xlim([0 16])
grid on 
fixticks
