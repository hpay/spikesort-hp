%% Generate list of directories to run kilosort on
birdID = 'HC05'; % AMB111
T = readtable('D:\hannah\Dropbox\alab\Analysis\RECORDING_DEPTH_CHICK.xlsx');
T = T(~T.exclude,:);
T = T(strcmp(T.bird, birdID),:)

% Data folder
dataDir = 'Z:\Hannah\Ephys\Project2';

%% set up code folders
scripts_dir = cd;
spikesort_hp_dir = fileparts(scripts_dir);
code_dir = fileparts(spikesort_hp_dir);
addpath(genpath(fullfile(spikesort_hp_dir,'src'))) % Add this second so any over-written functions are called e.g. get_session_waveforms
warning("off","parallel:gpu:device:DeviceDeprecated");

%% Add kilosort2 directory
addpath(genpath(fullfile(code_dir, 'kilosort-2.0')))

%% Run a single session using settings
%{
ops = [];
ops.chanMap = 'H6_shankA.mat';
ops = hp_config(ops);

rootDir = 'Z:\Hannah\ephys\project2\HC05_220819\raw_220819_125309'
saveDir = fullfile(rootDir,'kilo2.0'); mkdir(saveDir)
fprintf('\n Running Kilosort on directory %s \n', rootDir),
run_single_kilosort(rootDir,saveDir,ops)

% Save waveforms for celltype clustering
wvStruct = get_session_waveforms(rootDir, saveDir, 0);
save(fullfile(saveDir, 'waveformStruct.mat'), 'wvStruct')

% Load spikes into my dat structure in matlab
option_only_good= 1;
S = importKilo2(saveDir, ops.fs, option_only_good)

% Run phy to inspect it! (TODO)
% In Anaconda prompt:
% conda activate phy_env
% Navigate to folder with kilosort results (Type Z: to change to engram drive. then cd path. )
% phy template-gui params.py
%}

%% Run all sessions using settings in ops
for ii = 1:height(T)

    % Find the binary directory
    rootDir_temp = dir(fullfile(dataDir, T.filename{ii}, 'raw*'));
    rootDir = fullfile(rootDir_temp.folder, rootDir_temp.name);
        
    % Create save directory
    saveDir = fullfile(rootDir,'kilo2.0');  mkdir(saveDir);
         
    % set up options
    ops = [];
    ops.badChannels = [];
    ops.chanMap = [T.probe_chanmap{ii} '.mat']; % Make sure config folder is on the search path
    ops = hp_config(ops);    
    
    % Run and save results
    fprintf('\n Running Kilosort on directory %s \n', rootDir),
    run_single_kilosort(rootDir,saveDir, ops)
    
    % Save current options
    save(fullfile(saveDir,'ops.mat'),'ops')
    
    % Save waveforms for celltype clustering
    wvStruct = get_session_waveforms(rootDir, saveDir, 0);
    save(fullfile(saveDir, 'waveformStruct.mat'), 'wvStruct')
end

