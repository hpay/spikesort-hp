%% set up code folders
scripts_dir = cd;
spikesort_hp_dir = fileparts(scripts_dir);
code_dir = fileparts(spikesort_hp_dir);
addpath(genpath(fullfile(spikesort_hp_dir,'src'))) % Add this second so any over-written functions are called e.g. get_session_waveforms
warning("off","parallel:gpu:device:DeviceDeprecated");

%% Add ironclust directory
addpath(genpath(fullfile(code_dir, 'ironclust')))

%% Run a single session 
rootDir = 'Z:\Hannah\ephys\project2\HC05_220819\raw_220819_125309';
saveDir = fullfile(rootDir,'irc2'); mkdir(saveDir)
fprintf('\n Running Ironsort on directory %s \n', rootDir),
probeDir = which('H6_shankA.prb');

% Run ironclust with default settings
irc2(fullfile(rootDir,'amplifier.dat'), probeDir, saveDir)

% Export to Phy - DOESN'T WORK, use spikeinterface from here on! 
% Note: .mda files were introduced by mountainsort
% https://github.com/flatironinstitute/ironclust/issues/55
% irc2('export-phy', saveDir)

