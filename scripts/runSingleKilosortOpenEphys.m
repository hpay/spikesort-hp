function S = runSingleKilosortOpenEphys(session, ks_subdir, opts)
% RUNSINGLEKILOSORTOPENEPHYS  Spike sort one openEphys session with Kilosort 2 or 4, then waveforms, cell
% types (E / I / nonsomatic / noise, from the pooled whole-waveform model fit by fitCellTypeModel.m) and
% quality metrics (good / mua / noise). Returns the good units' spikes (importKilosort).
%
%   S = runSingleKilosortOpenEphys                    % defaults below
%   S = runSingleKilosortOpenEphys('MH02_20261111', 'kilosort2_output')
%   S = runSingleKilosortOpenEphys(session, ks_subdir, 'bad_channels_raw', [25 40], 'save_png', true)
%
% NOTHING is written to the server. Step 1 copies the session (including the raw continuous.dat,
% excluding video) to data_root, and every later step reads the local copy. Each step is skipped if its
% output already exists, so re-running just redoes the plots.
%
% ks_subdir: sort folder in the session (name contains kilosort2 or kilosort4; step 4 runs that Kilosort
%   if the folder has no sort yet). Current sorts:
%     kilosort4_mirrored_interleaved_nodrift   F_D_12shank_mirrored_interleaved.mat (current)
%     kilosort4_mirrored_aligned_nodrift       F_D_12shank_mirrored_aligned.mat
%     kilosort2_output                         F_D_12shank.mat
%     kilosort4 (1121), kilosort4_12shank_1200uM (1111)   original sorts with drift correction
% Channel map: an existing sort uses the map it was sorted with (findChanMap); chan_map_file is only for
%   new sorts. Maps are in src/kilosort_config_files, Cambridge NeuroTech F double-sided. Interleaved:
%   every face 200 um from the next (fronts at 400 um pitch, each back 200 um right of the front it pairs
%   with), so Kilosort treats each face as independent. Mirrored: back of A pairs with front of F, B with
%   E, C with D (the data correlate that way; the vendor sheet implies A with A). Aligned: each back at
%   the same x and y as its front. F_D_12shank.mat: backs at x 0-1000, fronts at 1200-2200.

arguments
    session = 'MH02_20261121'
    ks_subdir = 'kilosort4_mirrored_interleaved_nodrift'
    opts.chan_map_file = 'F_D_12shank_mirrored_interleaved.mat'   % for new sorts only (see above)
    opts.bad_channels_raw = 25      % raw openEphys channel numbers (CH#, 1-indexed) to exclude; see the channel noise plot (step 2)
    opts.bad_channels_auto = true   % also exclude the channels step 2 suggests
    opts.ks4_settings = struct('nblocks', 0)   % Kilosort 4 settings other than defaults (nblocks = 0: no drift correction); whitening_range is set in step 4
    opts.cell_type_model = 'cell_type_model_VW.mat'   % in src (or a full path), from fitCellTypeModel.m (VW = visual wulst); fit in step 6 if missing
    opts.overwrite_types = true     % re-apply the cell type model (step 6) even if gmm_result.mat exists
    % Step 7: which units are good (the rest of the E / I / nonsomatic units are mua)
    opts.min_type_p = 0.8           % cell type posterior below this: not good
    opts.min_spikes_good = 50       % this many spikes or fewer: not good (contamination: Kilosort 4's good label)
    opts.min_fit_p = 0              % E / I units whose fit to their E / I component has p below this: not good (0 = off)
    opts.save_png = false           % save figures as PNGs in <session>\figures (slow)
    opts.server_root = 'Z:\Meghana\Data\Ephys'   % READ ONLY
    opts.data_root = 'D:\dataMH'    % all outputs
end
server_dir = fullfile(opts.server_root, session);
local_dir = fullfile(opts.data_root, session);
copy_exceptions = {'\.mp4$'};   % don't copy video

%% Paths
code_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(code_dir, 'kilosort-2.0')))
addpath(genpath(fullfile(code_dir, 'src')))   % added last so src wins over older copies in kilosort-2.0 (e.g. getSpikeWaveform)
ks_dir = fullfile(local_dir, ks_subdir);
fig_dir = fullfile(local_dir, 'figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end
is_sorted = isfile(fullfile(ks_dir, 'spike_times.npy'));
chan_map_file = opts.chan_map_file;
if is_sorted
    chan_map_file = findChanMap(ks_dir);
    if isempty(chan_map_file); error('No channel map in src/kilosort_config_files matches %s', ks_dir); end
end
fprintf('%s %s: channel map %s\n', session, ks_subdir, chan_map_file)

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
[~, bad_suggested] = plotChannelNoise(raw.fbinary, raw.nCh, raw.fs, chan_map_file, raw.bit_volts(1), [], [], session);
if opts.save_png; exportgraphics(gcf, fullfile(fig_dir, 'channel_noise.png')); end
fprintf('Suggested bad channels (raw CH#): %s\n', mat2str(bad_suggested))
bad_channels = opts.bad_channels_raw;
if opts.bad_channels_auto; bad_channels = unique([bad_channels bad_suggested]); end
raw.chan_idx = setdiff(raw.ephys_idx, bad_channels); % bad channels are also left out of the waveforms (step 5)

%% 3. Threshold ADC channels (digital inputs) at 2.5 V and save rise/fall times (DIN.mat, same format as Intan)
dig_file = fullfile(local_dir, 'DIN.mat');
thresh_V = 2.5;
if ~exist(dig_file, 'file')
    thresh = round(thresh_V / raw.bit_volts(raw.adc_idx(1))); % in int16 counts
    DIN = importOpenEphysDigital(raw.fbinary, raw.nCh, raw.adc_idx, thresh, raw.fs, raw.channel_names(raw.adc_idx));
    save(dig_file, 'DIN')
end
% DIN = getfield(load(dig_file), 'DIN');
% plotDigitalEdges(DIN, 30); % opens on whole session; slider scrolls a 30 s window
% sgtitle([session ' digital inputs (>' num2str(thresh_V) ' V)'], 'Interpreter', 'none')

%% 4. Run Kilosort 2 or 4
if ~is_sorted && contains(ks_subdir, 'kilosort4')
    ks4_settings = opts.ks4_settings;
    ks4_settings.n_chan_bin = raw.nCh;
    ks4_settings.fs = raw.fs;
    m = load(chan_map_file);
    ks4_settings.whitening_range = min(32, min(accumarray(m.kcoords, 1))); % whiten within a face: fewest channels on any face (ignoring bad channels), max 32
    t = tic;
    runKilosort4(raw.fbinary, chan_map_file, ks_dir, bad_channels, ks4_settings);
    toc(t)
elseif ~is_sorted
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

%% 5. Extract waveforms for cell type clustering
if ~exist(fullfile(ks_dir, 'waveformStruct.mat'), 'file')
    only_good = 0; % all units; good is decided in step 7
    align_waveforms = true; % align each unit to its largest extremum (trough, or peak if bigger)
    wvStruct = getSessionWaveforms(local_dir, ks_dir, only_good, raw, align_waveforms);
    save(fullfile(ks_dir, 'waveformStruct.mat'), 'wvStruct', '-v7.3')
end

%% 6. Cell type of every unit (E / I / nonsomatic / noise) from the pooled whole-waveform model, and its
% posterior (type_p), saved in gmm_result.mat. Units without a usable waveform are noise. If the model file
% doesn't exist yet, it is fit first (fitCellTypeModel on its default sessions of this sort type).
cell_type_model = opts.cell_type_model;
model_file = cell_type_model; if ~isfile(model_file); model_file = which(cell_type_model); end   % full path, or a model in src
if isempty(model_file)
    model_file = fullfile(code_dir, 'src', cell_type_model);
    fprintf('No cell type model %s: fitting it\n', cell_type_model)
    fitCellTypeModel(ks_subdir, 'model_file', model_file);
end
model = load(model_file);
if ~exist(fullfile(ks_dir, 'gmm_result.mat'), 'file') || opts.overwrite_types
    [~, types, ~, type_ids, ~, res] = applyCellTypeModel(model.gm, {ks_dir}, 0, 0, 'waveform', model.basis);
    wvStruct = getfield(load(fullfile(ks_dir, 'waveformStruct.mat')), 'wvStruct');
    labels = repmat({'noise'}, numel(wvStruct.goodIDs), 1);       % same order as wvStruct (all units)
    type_p = NaN(numel(wvStruct.goodIDs), 1);                     % posterior of the unit's type
    fit_p = NaN(numel(wvStruct.goodIDs), 1);                      % E / I units: fit to their E / I component (cellTypesFromClusters)
    [found, loc] = ismember(wvStruct.goodIDs(:), type_ids); labels(found) = types(loc(found));
    type_p(found) = res.type_p(loc(found));
    if ~isempty(res.ei_fit_p); fit_p(found) = res.ei_fit_p(loc(found)); end
    id = wvStruct.goodIDs(:);
    save(fullfile(ks_dir, 'gmm_result.mat'), 'id', 'labels', 'type_p', 'fit_p', 'cell_type_model')
    wvStruct.typeLabels = labels;
    save(fullfile(ks_dir, 'waveformStruct.mat'), 'wvStruct', '-v7.3')
end

%% 7. Good / mua / noise labels from the cell types and the step 7 settings (applyQualityMetrics; skips
% phy-curated sessions), written to cluster_group.tsv and waveformStruct.mat; then the final cell type
% figure with those labels (orange x: cell type posterior below min_type_p; black +: E / I fit p below min_fit_p)
applyQualityMetrics(ks_dir, opts.min_type_p, opts.min_spikes_good, opts.min_fit_p);
if ~exist('res', 'var'); [~, ~, ~, ~, ~, res] = applyCellTypeModel(model.gm, {ks_dir}, 0, 0, 'waveform', model.basis); end
wvStruct = getfield(load(fullfile(ks_dir, 'waveformStruct.mat')), 'wvStruct');
[~, loc] = ismember(res.ids, wvStruct.goodIDs); res.unit_labels = wvStruct.goodLabels(loc);   % final labels
plotCellTypeModel(model.gm, res, 1, opts.min_type_p, opts.min_fit_p)
if opts.save_png; exportgraphics(gcf, fullfile(fig_dir, ['cell_types_' ks_subdir '.png'])); end

%% 8. Templates of the final good units at their positions on the probe
plotTemplatesOnProbe(ks_dir, [session ' ' ks_subdir ' (' chan_map_file ')']);
if opts.save_png; exportgraphics(gcf, fullfile(fig_dir, ['templates_on_probe_' ks_subdir '.png'])); end

%% 9. Load spikes into dat structure
option_only_good = 1;
S = importKilosort(ks_dir, raw.fs, option_only_good);

%% Run phy to inspect it!
% params.py points dat_path at the continuous.dat Kilosort read, with n_channels_dat = 136.
% In Anaconda prompt:
%   conda activate phy_env
%   D:
%   cd D:\dataMH\MH02_20261121\kilosort4_mirrored_interleaved_nodrift
%   phy template-gui params.py
end
