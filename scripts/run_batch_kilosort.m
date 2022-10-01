%% Generate list of directories to run kilosort on
birdID = 'HC05'; % AMB111
T = readtable('D:\hannah\Dropbox\alab\Analysis\RECORDING_DEPTH_CHICK.xlsx');
T = T(~T.exclude,:);
T = T(strcmp(T.bird, birdID),:)
fs = 3e4;

% Data folder
dataDir = 'Z:\Hannah\Ephys\Project2';

%% set up code folders
scripts_dir = cd;
spikesort_hp_dir = fileparts(scripts_dir);
code_dir = fileparts(spikesort_hp_dir);
addpath(genpath(fullfile(spikesort_hp_dir,'src'))) 
warning("off","parallel:gpu:device:DeviceDeprecated");

%% Add kilosort2 directory
addpath(genpath(fullfile(code_dir, 'kilosort-2.0')))

%% Run a single session using settings
% %{
ops = [];
% ops.chanMap = 'H6_shankA.mat';
ops.chanMap = 'H6.mat';
ops = hp_config(ops);

% rootDir = 'Z:\Hannah\ephys\project2\HC05_220819\raw_220819_125309'
rootDir = 'Z:\Hannah\ephys\project2\test64b\raw'
saveDir = fullfile(rootDir,'kilo2.0'); mkdir(saveDir)
fprintf('\n Running Kilosort on directory %s \n', rootDir)
t = tic
run_single_kilosort(rootDir,saveDir,ops)
toc(t)
% Save waveforms for celltype clustering
wvStruct = getSessionWaveforms(rootDir, saveDir, 0);
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
    wvStruct = getSessionWaveforms(rootDir, saveDir, 0);
    save(fullfile(saveDir, 'waveformStruct.mat'), 'wvStruct')
end


%% Collect # of good cells in each?
clear S
for ii = 1:height(T)


    % Find the binary directory
    rootDir_temp = dir(fullfile(dataDir, T.filename{ii}, 'raw*'));
    rootDir = fullfile(rootDir_temp.folder, rootDir_temp.name);
        
    % Create save directory
    saveDir = fullfile(rootDir,'kilo2.0');  
    
    option_only_good = 1;
    S{ii} =  importKilo2(saveDir, fs, option_only_good)
        
end
ncells = cellfun(@length,S)
days_since_surgery = days(T.date-datetime(2022,08,09))
% figure;
hold on
plot(days_since_surgery, ncells,'-ob')
xlabel('Days since surgery'); ylabel('# "good" cells in Kilosort2 (SC params)')
ylim([0 max(ncells)+2])
xlim([0 16])
grid on 
fixticks
