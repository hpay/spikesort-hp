%% Spike sort a single openEphys session with Kilosort 2.0
% Like runSingleKilosortExample.m, plus the waveform / GMM E-I / quality-metric
% steps from runBatchKilosort.m + analyzeBatchKilosort.m.
%
% NOTHING is written to the server. Step 1 copies the session (including the raw
% continuous.dat, excluding video) to local_dir, and every later step (ADC extraction,
% Kilosort, waveforms) reads the local copy. All outputs go to local_dir. Each step is
% skipped if its output already exists, so the whole script can be re-run to just redo
% the plots.

%% Settings
session = 'MH02_20261118';
% session = 'MH02_20261111';
server_dir = fullfile('Z:\Meghana\Data\Ephys', session);   % READ ONLY
local_dir = fullfile('D:\dataMH', session);                 % all outputs
chan_map_file = 'F_D_12shank.mat'; % Cambridge NeuroTech F double-sided, back treated as offset in x (src/kilosort_config_files)
bad_channels_raw = [];  % raw openEphys channelanumbers (CH#, 1-indexed) to exclude. See channel noise plot (step 3)
bad_channels_auto = true;
ks_subdir = 'kilosort2_output'; % sort folder to process (name contains kilosort2 or kilosort4): Kilosort 2.0 runs here in step 4; a KS4 sort must already exist, e.g. 'kilosort4_12shank_1200uM'
% ks_subdir = 'kilosort4_12shank_1200uM'; % sort folder to process (name contains kilosort2 or kilosort4): Kilosort 2.0 runs here in step 4; a KS4 sort must already exist, e.g. 'kilosort4_12shank_1200uM'
overwrite_gmm = 1;      % set 1 to retrain GMM
gmm_label_by = 'width'; % name GMM clusters E/I by spike width (narrow = I); 'rate' = lower rate is E

% Don't copy video
copy_exceptions = {'\.mp4$'};

%% Paths
code_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(code_dir, 'src')))
addpath(genpath(fullfile(code_dir, 'kilosort-2.0')))
ks_dir_fun = @(root_dir) fullfile(root_dir, ks_subdir);
ks_dir = ks_dir_fun(local_dir);
gmm_file = fullfile(ks_dir, 'gmm_model.mat'); % GMM trained on this sort only
fig_dir = fullfile(local_dir, 'figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

%% 1. Copy the session (except video) to local_dir, then read raw data locally
% Files already copied are skipped, so re-running only reads the server file listing.
f = dir(fullfile(server_dir, '**', '*'));
f = f(~[f.isdir] & cellfun(@isempty, regexp({f.name}, strjoin(copy_exceptions, '|'))));
need_GB = sum([f(~isfile(fullfile(strrep({f.folder}, server_dir, local_dir), {f.name}))).bytes])/1e9; % files still to copy
free_GB = java.io.File(local_dir(1:3)).getUsableSpace()/1e9;
if free_GB < need_GB; error('Need %.1f GB free on %s, have %.1f GB', need_GB, local_dir(1:3), free_GB); end
runmode = 2; % never delete, never overwrite newer files
dirbackup(server_dir, local_dir, runmode, copy_exceptions)

local_oebin = dir(fullfile(local_dir, '**', 'structure.oebin'));
if numel(local_oebin) ~= 1; error('Expected one local structure.oebin, found %d', numel(local_oebin)); end
raw = readOpenEphysInfo(local_oebin.folder);
raw.chan_idx = raw.ephys_idx; % CH1-128 = amplifier channels; ADC1-8 after that
fprintf('%s: %d channels (%d ephys, %d ADC), %.1f min at %g Hz\nReading raw data from %s\n', session, raw.nCh, ...
    numel(raw.ephys_idx), numel(raw.adc_idx), raw.nSamp/raw.fs/60, raw.fs, raw.fbinary)

%% 2. Threshold ADC channels (digital inputs) at 2.5 V and save rise/fall times (DIN.mat, same format as Intan)
dig_file = fullfile(local_dir, 'DIN.mat');
thresh_V = 2.5;
if ~exist(dig_file, 'file')
    thresh = round(thresh_V / raw.bit_volts(raw.adc_idx(1))); % in int16 counts
    DIN = importOpenEphysDigital(raw.fbinary, raw.nCh, raw.adc_idx, thresh, raw.fs, raw.channel_names(raw.adc_idx));
    save(dig_file, 'DIN')
end
DIN = getfield(load(dig_file), 'DIN');
plotDigitalEdges(DIN, 30); % opens on whole session; slider scrolls a 30 s window
sgtitle([session ' digital inputs (>' num2str(thresh_V) ' V)'], 'Interpreter', 'none')
exportgraphics(gcf, fullfile(fig_dir, 'digital_inputs.png'))

%% 3. Channel noise QC (pick bad_channels_raw from this)
[noise_uV, bad_suggested] = plotChannelNoise(raw.fbinary, raw.nCh, raw.fs, chan_map_file, raw.bit_volts(1));
sgtitle(session, 'Interpreter', 'none')
exportgraphics(gcf, fullfile(fig_dir, 'channel_noise.png'))
fprintf('Suggested bad channels (raw CH#): %s\n', mat2str(bad_suggested))

%% 4. Run Kilosort 2.0
if isempty(dir(fullfile(ks_dir, 'spike_times.npy')))
    if ~contains(ks_subdir, 'kilosort2'); error('No sort in %s (only Kilosort 2.0 runs here)', ks_dir); end
    if ~exist(ks_dir, 'dir'); mkdir(ks_dir); end
    ops = [];
    ops.chanMap = chan_map_file;
    ops = hp_config(ops);
    ops.NchanTOT = raw.nCh;     % 136: KS reads only the chanMap channels
    ops.fs = raw.fs;
    ops.fbinary = raw.fbinary;  % local copy
    m = load(chan_map_file);
    if bad_channels_auto; bad_channels = [bad_channels_raw bad_suggested]; else
        bad_channels = bad_channels_raw; end
    ops.badChannels = find(ismember(m.chanMap, bad_channels)); % KS indexes bad channels by position in chanMap
    if ~exist(fileparts(ops.fproc), 'dir'); mkdir(fileparts(ops.fproc)); end

    fprintf('\n Running Kilosort on %s \n', ops.fbinary)
    t = tic;
    runSingleKilosort(local_dir, ks_dir, ops);
    toc(t)
    delete(ops.fproc) % ~40 GB whitened temp file
else
    fprintf('Already sorted %s\n', ks_dir)
end

%% 5. Waveforms for cell type clustering
if ~exist(fullfile(ks_dir, 'waveformStruct.mat'), 'file')
    only_good = 0; % Process all for now, select good later
    align_waveforms = true; % align each unit to its largest extremum (trough, or peak if bigger)
    wvStruct = getSessionWaveforms(local_dir, ks_dir, only_good, raw, align_waveforms);
    save(fullfile(ks_dir, 'waveformStruct.mat'), 'wvStruct', '-v7.3')
end

%% 6. Train a GMM on this session (Kilosort's original "good" units), then apply to all units
if ~exist(gmm_file, 'file') || overwrite_gmm % Train model
    option_only_good = 2; % Kilosort's labels, not those from step 7 (which uses the GMM output)
    n_clusters = 2;
    gm = GMM_make({ks_dir}, option_only_good, n_clusters, gmm_label_by);
    exportgraphics(gcf, fullfile(fig_dir, ['gmm_train_' ks_subdir '.png']))
    save(gmm_file, 'gm', 'session')
    delete(fullfile(ks_dir, 'gmm_result.mat')) % re-apply the new model below
end
if ~exist(fullfile(ks_dir, 'gmm_result.mat'), 'file') % Apply model
    gm = getfield(load(gmm_file), 'gm');
    [id, labels, stats] = GMM_apply(gm, ks_dir, 0, 1);
    exportgraphics(gcf, fullfile(fig_dir, ['gmm_apply_' ks_subdir '.png']))
    save(fullfile(ks_dir, 'gmm_result.mat'), 'id', 'labels', 'stats')
    wvStruct = getfield(load(fullfile(ks_dir, 'waveformStruct.mat')), 'wvStruct');
    wvStruct.typeLabels = labels;
    save(fullfile(ks_dir, 'waveformStruct.mat'), 'wvStruct')
end

%% 7. Basic sorting quality metrics: updates good/mua/noise labels (skips phy-curated sessions)
applyQualityMetrics(ks_dir);

%% 8. Plot E/I clusters
T = table({session}, 'VariableNames', {'filename'});
plotEIClusters(T, fileparts(local_dir), ks_dir_fun, false, 1, 2)
title([session ' ' ks_subdir], 'Interpreter', 'none')
exportgraphics(gcf, fullfile(fig_dir, ['EI_clusters_' ks_subdir '.png']))

%% 9. Load spikes into dat structure
option_only_good = 1;
S = importKilosort(ks_dir, raw.fs, option_only_good);
disp(S)

%% Run phy to inspect it!
% params.py points dat_path at the continuous.dat Kilosort read, with n_channels_dat = 136.
% In Anaconda prompt:
%   conda activate phy_env
%   D:
%   cd D:\dataMH\MH02_20261118\kilosort2_output
%   phy template-gui params.py
