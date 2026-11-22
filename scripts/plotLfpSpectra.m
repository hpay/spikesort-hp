%% Multitaper power spectrum of each channel vs depth below the brain surface, pooled across
% sessions, one panel per shank.
% Spectra: Chronux mtspectrumc on 300 x 2 s segments spread through each session, downsampled
% to 1 kHz, tapers [TW K] = [2 3] (1 Hz half-bandwidth), averaged over segments (see channelSpectra).
% Short windows keep oscillations whose frequency wanders from being split across narrow bins.
% Each session's spectra are saved in its folder (lfp_spectra.mat); the figure goes to data_dir\figures.

sessions = {'MH02_20261111', 'MH02_20261118'};
data_dir = 'D:\dataMH';
chan_map_file = 'F_D_12shank.mat';
opts = struct('n_seg', 300, 'seg_dur', 2, 'tapers', [2 3]);

code_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(code_dir, 'src')))
fig_dir = fullfile(data_dir, 'figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

%% Compute (or load) each session's spectra, in dB relative to the shank's median channel at each frequency
m = load(chan_map_file);
power_dB = []; depth = []; shank = [];
for ii = 1:numel(sessions)
    local_dir = fullfile(data_dir, sessions{ii});
    spectra_file = fullfile(local_dir, 'lfp_spectra.mat');
    if ~exist(spectra_file, 'file')
        raw = readOpenEphysInfo(getfield(dir(fullfile(local_dir, '**', 'structure.oebin')), 'folder'));
        t = tic;
        [S, f] = channelSpectra(raw.fbinary, raw.nCh, raw.fs, m.chanMap, raw.bit_volts(1), opts);
        fprintf('%s: spectra computed in %.0f s\n', sessions{ii}, toc(t))
        chan_map = m.chanMap; % columns of S are these raw CH#
        save(spectra_file, 'S', 'f', 'chan_map', 'opts')
    end
    load(spectra_file, 'S', 'f', 'chan_map')
    [depth_ii, shank_ii] = channelDepths(sessions{ii}, chan_map, data_dir, chan_map_file);
    p = 10*log10(S);
    for k = 1:12 % normalize within each shank and session
        p(:, shank_ii == k) = p(:, shank_ii == k) - median(p(:, shank_ii == k), 2);
    end
    power_dB = [power_dB, p]; %#ok<AGROW>
    depth = [depth; depth_ii]; %#ok<AGROW>
    shank = [shank; shank_ii]; %#ok<AGROW>
end

%% Plot: power vs frequency and depth, 3 x 4 shanks (depths with no channel left blank)
row_pitch = median(diff(unique(m.ycoords)));
depth_grid = min(depth) : row_pitch : max(depth);
c_lim = [-1 1] * prctile(abs(power_dB(:)), 99);
figure;
for k = 1:12
    img = nan(numel(f), numel(depth_grid));
    on_shank = find(shank == k);
    [~, row] = min(abs(depth(on_shank) - depth_grid), [], 2);
    img(:, row) = power_dB(:, on_shank);
    subplot(3, 4, k)
    imagesc(f, depth_grid, img', 'AlphaData', ~isnan(img'))
    clim(c_lim); colormap(gca, redblueColormap()); set(gca, 'TickDir', 'out', 'Box', 'off')
    title(sprintf('shank %d', k))
    if mod(k-1, 4) == 0; ylabel('depth below surface (\mum)'); end
    if k > 8; xlabel('frequency (Hz)'); end
end
c = colorbar('Position', [0.93 0.11 0.012 0.8]); c.Label.String = 'power re shank median (dB)';
sgtitle({sprintf('LFP power by channel: %s', strjoin(sessions, ', ')), ...
    're median channel of each shank and session, at each frequency'}, 'Interpreter', 'none')
drawnow
exportgraphics(gcf, fullfile(fig_dir, 'lfp_spectra_by_depth.png'))

function cmap = redblueColormap()
% blue - white - red, 256 levels
x = linspace(0, 1, 128)';
cmap = [[x x ones(128,1)]; [ones(128,1) flipud(x) flipud(x)]];
end
