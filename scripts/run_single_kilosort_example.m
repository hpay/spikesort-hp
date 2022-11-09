%% Run a single session 

% Select file and probe
raw_dir = 'Z:\Hannah\ephys\project2\HC09_221104\raw_221104_121402';
ops = []; ops.chanMap = 'H6.mat';

% raw_dir = 'Z:\Hannah\ephys\project2\HC05_220819\raw_220819_125309'
% ops = []; ops.chanMap = 'H6_shankA.mat';

% Run sorting
ops = hp_config(ops);
save_dir = fullfile(fileparts(raw_dir,'kilosort2_output')); 
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
S = importKilosort(save_dir, ops.fs, option_only_good);
disp(S)

% Run phy to inspect it! (TODO)
% In Anaconda prompt:
%   conda activate phy_env
% Navigate to folder with kilosort results (Type Z: to change to engram drive. then cd path. )
%   phy template-gui params.py
%}