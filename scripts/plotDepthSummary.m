%% Summary of sorted units as a function of depth below the brain surface, pooled across sessions
% Depth = probe position in the brain (DataNotebook-MH.xlsx, Ephys sheet) - channel y (height above the tip)
% Uses KS "good" units (after applyQualityMetrics) with GMM label E or I.

sessions = {'MH02_20261111', 'MH02_20261118'};
ks_subdir = 'kilosort2_output';
chan_map_file = 'F_D_12shank.mat';
data_dir = 'D:\dataMH';

code_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(code_dir, 'src')))
fig_dir = fullfile(data_dir, 'figures');
if ~exist(fig_dir, 'dir'); mkdir(fig_dir); end

%% Load units: depth of the channel with each unit's largest mean waveform, and cell type
depth = []; cell_type = {}; ks_label = {};
for ii = 1:numel(sessions)
    wv = getfield(load(fullfile(data_dir, sessions{ii}, ks_subdir, 'waveformStruct.mat')), 'wvStruct');
    % max_site: raw CH# with the largest peak-to-peak mean waveform (NaN if no spikes)
    [unit_depth, ~, tip_depth] = channelDepths(sessions{ii}, wv.max_site, data_dir, chan_map_file);
    labels = wv.typeLabels(:); good_mua = wv.goodLabels(:);
    keep = ~isnan(unit_depth) & ismember(good_mua, {'good', 'mua'}) & ismember(labels, {'E', 'I'});
    depth = [depth; unit_depth(keep)]; %#ok<AGROW> depth below the brain surface
    cell_type = [cell_type; labels(keep)]; %#ok<AGROW>
    ks_label = [ks_label; good_mua(keep)]; %#ok<AGROW>
    fprintf('%s: tip %g um below the surface; good E %d, I %d; MUA E %d, I %d\n', sessions{ii}, tip_depth, ...
        nnz(keep & strcmp(labels, 'E') & strcmp(good_mua, 'good')), nnz(keep & strcmp(labels, 'I') & strcmp(good_mua, 'good')), ...
        nnz(keep & strcmp(labels, 'E') & strcmp(good_mua, 'mua')), nnz(keep & strcmp(labels, 'I') & strcmp(good_mua, 'mua')))
end

%% 1. Fraction of I cells, I/(E+I), vs depth, with 95% binomial (Clopper-Pearson) CIs: good, and good + MUA
depth_bin = 15; % um; 15 = one row of sites
bin_edges = min(depth) - depth_bin/2 : depth_bin : max(depth) + depth_bin; % last edge past the deepest unit
bin_centers = bin_edges(1:end-1) + depth_bin/2;
bin = discretize(depth, bin_edges);
unit_sets = {{'good'}, {'good', 'mua'}};
set_names = {'good', 'good + MUA'};
colors = {'k', [0.85 0.2 0.2]};
y_offset = [-2 2]; % um, so the two sets don't overlap

figure; hold on
h = gobjects(2,1); total = zeros(2,2);
for s = 1:2
    in = ismember(ks_label, unit_sets{s});
    n_I = accumarray(bin(in), strcmp(cell_type(in), 'I'), [numel(bin_centers) 1]);
    n_EI = accumarray(bin(in), 1, [numel(bin_centers) 1]);
    has_cells = n_EI > 0;
    [frac_I, ci] = binofit(n_I(has_cells), n_EI(has_cells));
    y = bin_centers(has_cells)' + y_offset(s);
    h(s) = errorbar(frac_I, y, frac_I - ci(:,1), ci(:,2) - frac_I, 'horizontal', 'o', 'Color', colors{s}, ...
        'MarkerFaceColor', colors{s}, 'CapSize', 3);
    text(ci(:,2) + 0.02, y, compose('%d/%d', n_I(has_cells), n_EI(has_cells)), 'FontSize', 7, 'Color', colors{s})
    total(s,:) = [sum(n_I) sum(n_EI)];
end
set(gca, 'YDir', 'reverse', 'Box', 'off', 'TickDir', 'out')
xlim([0 1]); ylim(bin_edges([1 end]))
xlabel('fraction I cells, I/(E+I)'); ylabel('depth below brain surface (\mum)')
legend(h, compose('%s (%d I / %d E+I)', string(set_names'), total(:,1), total(:,2)), 'Location', 'east', 'Box', 'off')
title({strjoin(sessions, ', '), 'all shanks'}, 'Interpreter', 'none')
drawnow
exportgraphics(gcf, fullfile(fig_dir, ['depth_fraction_I_' ks_subdir '.png']))
