%% Generate list of directories to run kilosort on
% allDir = dir('D:\LMN164*');
birdID = 'HC05'; % AMB111
T = readtable('D:\hannah\Dropbox\alab\Analysis\RECORDING_DEPTH_CHICK.xlsx');
T = T(strcmp(T.bird, birdID),:)

% Data folder
dataDir = 'Z:\Hannah\Ephys\Project2';

%% set up code folders
scripts_dir = cd;
spikesort_hp_dir = fileparts(scripts_dir);
code_dir = fileparts(spikesort_hp_dir);
addpath(genpath(fullfile(code_dir, 'kilosort-2.0')))
addpath(genpath(fullfile(spikesort_hp_dir,'src'))) % Add this second so any over-written functions are called e.g. get_session_waveforms

%% set up channel map
ii = 1;
% chanMapFile = fullfile(spikesort_hp_dir, 'src\kilosort_config_files\', [T.kilosort_chanmap{ii} '.mat']);
chanMapFile = [T.kilosort_chanmap{ii} '.mat'];
badChannels = [];

%% set up options
ops = [];
ops.badChannels = [];
ops.chanMap = chanMapFile; % Chanmap
ops = hp_config(ops);

%% Run a single session using settings
rootDir = 'Z:\Hannah\ephys\project2\HC05_220819\raw_220819_125309'
saveDir = fullfile(rootDir,'kilo2.0'); mkdir(saveDir)
fprintf('\n Running Kilosort on directorytmp %s \n', rootDir),
run_single_kilosort(rootDir,saveDir,ops)

% Load spikes into my dat structure in matlab
option_only_good= 1;
S = importKilo2(saveDir, ops.fs, option_only_good)


% Run phy to inspect it! (TODO)
system(['phy template-gui ' fullfile(saveDir,'params.py')])


%% Run all sessions using settings in ops
for ii = height(T)
    
    % Find the binary directory
    rootDir_temp = dir(fullfile(dataDir, T.filename{ii}, 'raw*'));
    rootDir = fullfile(rootDir_temp.folder, rootDir_temp.name);
        
    % Create save directory
    saveDir = fullfile(rootDir,'kilo2.0'); % path to temporary binary file (same size as data, ideally on fast SSD)
    if ~exist(saveDir,'dir')
        mkdir(saveDir)
    end
            
    % Run and save results
    fprintf('\n Running Kilosort on directorytmp %s \n', rootDir),
    run_single_kilosort(rootDir,saveDir, ops)
    
    % Save current options
    save(fullfile(saveDir,'ops.mat'),'ops')
    
    % Save waveforms for celltype clustering
    wvStruct = get_session_waveforms(rootDir, saveDir, 0);
    fname = fullfile(saveDir, 'waveformStruct.mat');
    fprintf('Saving mean waveform data \n')
    save(fname, 'wvStruct')
end

