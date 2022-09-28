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

%% Add ironclust directory
addpath(genpath(fullfile(code_dir, 'ironclust')))

%% Run a single session using settings
% %{

rootDir = 'Z:\Hannah\ephys\project2\HC05_220819\raw_220819_125309';
saveDir = fullfile(rootDir,'irc2'); mkdir(saveDir)
fprintf('\n Running Ironsort on directory %s \n', rootDir),
probeDir = which('H6_shankA.prb');

% Run ironclust with default settings
irc2(fullfile(rootDir,'amplifier.dat'), probeDir, saveDir)

% Export to Phy - DOESN'T WORK, use spikeinterface from here on! .mda files
% were introduced by mountainsort
% https://github.com/flatironinstitute/ironclust/issues/55
% irc2('export-phy', saveDir)


%% Run all sessions 
for ii = 1:height(T)

    % Find the binary directory
    rootDir_temp = dir(fullfile(dataDir, T.filename{ii}, 'raw*'));
    rootDir = fullfile(rootDir_temp.folder, rootDir_temp.name);
        
    % Create save directory
    saveDir = fullfile(rootDir,'irc2'); mkdir(saveDir); 
         
    % Get probe directory
    probeDir = which([T.probe_chanmap{ii} '.prb']);

    % Run and save results
    fprintf('\n Running Ironclust on directory %s \n', rootDir),
    irc2(fullfile(rootDir,'amplifier.dat'), probeDir, saveDir)

end

