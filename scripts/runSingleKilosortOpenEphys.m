%% Spike sort a single openEphys session with Kilosort 2 or 4
% Like runSingleKilosortExample.m, plus waveforms, cell types (E / I / nonsomatic / noise, from the pooled
% whole-waveform model fit by fitCellTypeModel.m) and quality metrics (good / mua / noise).
%
% NOTHING is written to the server. Step 1 copies the session (including the raw
% continuous.dat, excluding video) to local_dir, and every later step (ADC extraction,
% Kilosort, waveforms) reads the local copy. All outputs go to local_dir. Each step is
% skipped if its output already exists, so the whole script can be re-run to just redo
% the plots.

%% Settings
% session = 'MH02_20261118';
% session = 'MH02_20261111';
session = 'MH02_20261121';
if exist('session_override', 'var'); session = session_override; end % lets a batch runner pick the session
server_dir = fullfile('Z:\Meghana\Data\Ephys', session);   % READ ONLY
local_dir = fullfile('D:\dataMH', session);                 % all outputs
% Channel maps (src/kilosort_config_files), Cambridge NeuroTech F double-sided. Interleaved: every face 200 um
% from the next (fronts at 400 um pitch, each back 200 um right of the front it pairs with), so Kilosort treats
% each face as independent. Mirrored: back of A pairs with front of F, B with E, C with D (the data correlate
% that way; the vendor sheet implies A with A).
% chan_map_file = 'F_D_12shank_mirrored_interleaved.mat';
% chan_map_file = 'F_D_12shank_mirrored_aligned.mat'; % mirrored, each back at the same x (and y) as its front
% chan_map_file = 'F_D_12shank_interleaved.mat'; % A_back with A_front. Looks wrong
chan_map_file = 'F_D_12shank.mat';             % original: backs at x 0-1000, fronts at 1200-2200 (used for kilosort2_output, kilosort4_nodrift)
if exist('chan_map_override', 'var'); chan_map_file = chan_map_override; end
bad_channels_raw = 25;  % raw openEphys channel numbers (CH#, 1-indexed) to exclude. See channel noise plot (step 3)
bad_channels_auto = true;
% ks_subdir = 'kilosort2_output';  % sorted with F_D_12shank.mat
% ks_subdir = 'kilosort2_mirrored'; % sort folder to process (name contains kilosort2 or kilosort4); step 4 runs the matching Kilosort if the folder has no sort yet
% ks_subdir = 'kilosort4_mirrored_interleaved_nodrift'; % with F_D_12shank_mirrored_interleaved.mat
% ks_subdir = 'kilosort4_mirrored_aligned_nodrift';     % with F_D_12shank_mirrored_aligned.mat
% ks_subdir = 'kilosort4_nodrift'; % sorted with F_D_12shank.mat
ks_subdir = 'kilosort4'; % sorted with F_D_12shank.mat - original, with drift correction
if exist('ks_subdir_override', 'var'); ks_subdir = ks_subdir_override; end
ks4_settings = struct('nblocks', 0); % Kilosort 4 settings other than defaults (nblocks = 0: no drift correction); whitening_range is set in step 4
cell_type_model = 'cell_type_model_VW.mat'; % cell type model (E / I / nonsomatic / noise) for the recorded brain region, in src, from fitCellTypeModel.m (VW = visual wulst)
overwrite_types = 1;    % set 1 to re-apply the cell type model (step 6) even if gmm_result.mat exists
save_png = false;       % save figures as PNGs in local_dir\figures (slow)
zero_index_figs = false; % label channels 0-indexed in the channel noise figure (probe maker / Intan numbering)

% Don't copy video
copy_exceptions = {'\.mp4$'};

%% Paths
code_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(code_dir, 'kilosort-2.0')))
addpath(genpath(fullfile(code_dir, 'src')))   % added last so src wins over older copies in kilosort-2.0 (e.g. getSpikeWaveform)
ks_dir_fun = @(root_dir) fullfile(root_dir, ks_subdir);
ks_dir = ks_dir_fun(local_dir);
fig_dir = fullfile(local_dir, 'figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

%% 1. Copy the session (except video) to local_dir, then read raw data locally
% Re-running only reads the server file listing; dirbackup runs only if a file is missing locally.
f = dir(fullfile(server_dir, '**', '*'));
f = f(~[f.isdir] & cellfun(@isempty, regexp({f.name}, strjoin(copy_exceptions, '|'))));
missing = ~isfile(fullfile(strrep({f.folder}, server_dir, local_dir), {f.name}));
if any(missing)
    need_GB = sum([f(missing).bytes])/1e9;
    free_GB = java.io.File(local_dir(1:3)).getUsableSpace()/1e9;
    if free_GB < need_GB; error('Need %.1f GB free on %s, have %.1f GB', need_GB, local_dir(1:3), free_GB); end
    runmode = 2; % never delete, never overwrite newer files
    dirbackup(server_dir, local_dir, runmode, copy_exceptions)
end

local_oebin = dir(fullfile(local_dir, '**', 'structure.oebin'));
if numel(local_oebin) ~= 1; error('Expected one local structure.oebin, found %d', numel(local_oebin)); end
raw = readOpenEphysInfo(local_oebin.folder);
raw.chan_idx = raw.ephys_idx; % CH1-128 = amplifier channels; ADC1-8 after that
fprintf('%s: %d channels (%d ephys, %d ADC), %.1f min at %g Hz\nReading raw data from %s\n', session, raw.nCh, ...
    numel(raw.ephys_idx), numel(raw.adc_idx), raw.nSamp/raw.fs/60, raw.fs, raw.fbinary)

%% 2. Channel noise QC (pick bad_channels_raw from this)
[noise_uV, bad_suggested] = plotChannelNoise(raw.fbinary, raw.nCh, raw.fs, chan_map_file, raw.bit_volts(1), [], [], session, zero_index_figs);
if save_png; exportgraphics(gcf, fullfile(fig_dir, 'channel_noise.png')); end
fprintf('Suggested bad channels (raw CH#): %s\n', mat2str(bad_suggested))
if bad_channels_auto; bad_channels = unique([bad_channels_raw bad_suggested]); else
    bad_channels = bad_channels_raw; end
raw.chan_idx = setdiff(raw.ephys_idx, bad_channels); % bad channels are also left out of the waveforms (step 5)

%% 3. Threshold ADC channels (digital inputs) at 2.5 V and save rise/fall times (DIN.mat, same format as Intan)
dig_file = fullfile(local_dir, 'DIN.mat');
thresh_V = 2.5;
if ~exist(dig_file, 'file')
    thresh = round(thresh_V / raw.bit_volts(raw.adc_idx(1))); % in int16 counts
    DIN = importOpenEphysDigital(raw.fbinary, raw.nCh, raw.adc_idx, thresh, raw.fs, raw.channel_names(raw.adc_idx));
    save(dig_file, 'DIN')
end
DIN = getfield(load(dig_file), 'DIN');
% plotDigitalEdges(DIN, 30); % opens on whole session; slider scrolls a 30 s window
% sgtitle([session ' digital inputs (>' num2str(thresh_V) ' V)'], 'Interpreter', 'none')
% if save_png; exportgraphics(gcf, fullfile(fig_dir, 'digital_inputs.png')); end

%% 4. Run Kilosort 2 or 4
if isempty(dir(fullfile(ks_dir, 'spike_times.npy'))) && contains(ks_subdir, 'kilosort4')
    ks4_settings.n_chan_bin = raw.nCh; 
    ks4_settings.fs = raw.fs;
    m = load(chan_map_file); 
    ks4_settings.whitening_range = min(32, min(accumarray(m.kcoords, 1))); % whiten within a face: fewest channels on any face (ignoring bad channels), max 32
    t = tic;
    runKilosort4(raw.fbinary, chan_map_file, ks_dir, bad_channels, ks4_settings);
    toc(t)
elseif isempty(dir(fullfile(ks_dir, 'spike_times.npy')))
    if ~exist(ks_dir, 'dir'); mkdir(ks_dir); end
    ops = [];
    ops.chanMap = chan_map_file;
    ops = hp_config(ops);
    ops.NchanTOT = raw.nCh;     % 136: KS reads only the chanMap channels
    ops.fs = raw.fs;
    ops.fbinary = raw.fbinary;  % local copy
    m = load(chan_map_file);
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

%% 6. Cell type of every unit (E / I / nonsomatic / noise) from the pooled whole-waveform model
% (fit across sessions by fitCellTypeModel.m). Units without a usable waveform are noise.
if ~exist(fullfile(ks_dir, 'gmm_result.mat'), 'file') || overwrite_types
    model = load(which(cell_type_model));
    [~, types, ~, type_ids] = GMM_apply(model.gm, {ks_dir}, 0, 1, 'waveform', model.basis);
    if save_png; exportgraphics(gcf, fullfile(fig_dir, ['cell_types_' ks_subdir '.png'])); end
    wvStruct = getfield(load(fullfile(ks_dir, 'waveformStruct.mat')), 'wvStruct');
    labels = repmat({'noise'}, numel(wvStruct.goodIDs), 1);       % same order as wvStruct (all units)
    [found, loc] = ismember(wvStruct.goodIDs(:), type_ids); labels(found) = types(loc(found));
    id = wvStruct.goodIDs(:);
    save(fullfile(ks_dir, 'gmm_result.mat'), 'id', 'labels', 'cell_type_model')
    wvStruct.typeLabels = labels;
    save(fullfile(ks_dir, 'waveformStruct.mat'), 'wvStruct', '-v7.3')
end

%% 7. Basic sorting quality metrics: updates good/mua/noise labels (skips phy-curated sessions)
applyQualityMetrics(ks_dir);

%% 8. Templates of the final good units at their positions on the probe
% (cell types are plotted in step 6)
plotTemplatesOnProbe(ks_dir, [session ' ' ks_subdir ' (' chan_map_file ')']);
if save_png; exportgraphics(gcf, fullfile(fig_dir, ['templates_on_probe_' ks_subdir '.png'])); end

%% 9. Load spikes into dat structure
option_only_good = 1;
S = importKilosort(ks_dir, raw.fs, option_only_good);
disp(S)

%% Run phy to inspect it!
% params.py points dat_path at the continuous.dat Kilosort read, with n_channels_dat = 136.
% In Anaconda prompt:
%   conda activate phy_env
%   D:
%   cd D:\dataMH\MH02_20261121\kilosort4_nodrift
%   phy template-gui params.py
