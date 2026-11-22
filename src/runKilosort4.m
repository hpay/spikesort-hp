function runKilosort4(data_file, chan_map_file, results_dir, bad_channels, ks4_settings)
% RUNKILOSORT4  Run Kilosort 4 (Python) from MATLAB, using the kilosort4 conda environment
% and the Kilosort source in spikesort-hp/kilosort4. Blocks until the sort finishes; Kilosort's
% console output is shown in the command window and its log is saved in results_dir.
%
%   data_file: raw int16 binary; chan_map_file: Kilosort2 .mat channel map (e.g. F_D_12shank.mat)
%   bad_channels: raw CH# (1-indexed) to exclude
%   ks4_settings: struct of Kilosort 4 settings (must include n_chan_bin and fs), e.g.
%       struct('n_chan_bin', 136, 'fs', 30000, 'nblocks', 0)

python = 'C:\Users\Hannah\.conda\envs\kilosort4\python.exe';
code_dir = fileparts(mfilename('fullpath'));

cfg = struct('data_file', data_file, 'probe_file', which(chan_map_file), 'results_dir', results_dir, ...
    'bad_channels', bad_channels(:)', 'ks4_settings', ks4_settings);
if ~exist(results_dir, 'dir'); mkdir(results_dir); end
cfg_file = fullfile(results_dir, 'ks4_run_settings.json');
writelines(jsonencode(cfg, 'PrettyPrint', true), cfg_file)

status = system(sprintf('"%s" "%s" "%s"', python, fullfile(code_dir, 'runKilosort4.py'), cfg_file), '-echo');
if status ~= 0; error('Kilosort 4 failed (exit status %d); see %s', status, fullfile(results_dir, 'kilosort4.log')); end
end
