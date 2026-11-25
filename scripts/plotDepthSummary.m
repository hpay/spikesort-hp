%% Good units of each cell type (E / I / nonsomatic) vs depth below the brain surface, by shank, pooled across sessions
% Depth = probe tip position in the brain (DataNotebook-MH.xlsx, Ephys sheet) - height above the tip of
% each unit's peak channel (wvStruct.max_site). Shanks are the physical shanks: the two faces of one shank
% pair in the mirrored map, named by the front face (its back is the mirrored partner, e.g. front A with
% back F). Cell types from step 6 and good labels from step 7 of runSingleKilosortOpenEphys.
% 1. Each good unit at its shank and depth, stacked like a histogram
% 2. Fraction of each type per depth bin, with 95% Clopper-Pearson CIs: each shank, then all shanks pooled
% Also prints each type's fraction per shank and pooled (all depths).

sessions = {'MH02_20261111', 'MH02_20261118', 'MH02_20261121'};
ks_subdir = 'kilosort4_mirrored_interleaved_nodrift';
data_dir = 'D:\dataMH';
depth_bin = 50;          % um
pair_pitch = 400;        % um: one physical shank per 400 um of map x (interleaved and aligned maps)
types = {'E', 'I', 'nonsomatic'};
type_color = {[0.85 0.1 0.1], [0.1 0.3 0.85], [0.1 0.6 0.2]};
font_size = 10;

code_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(code_dir, 'src')))

%% Load good units: depth and shank of the peak channel, cell type
depth = []; shank = []; cell_type = {};
for ii = 1:numel(sessions)
    ks_dir = fullfile(data_dir, sessions{ii}, ks_subdir);
    wv = getfield(load(fullfile(ks_dir, 'waveformStruct.mat')), 'wvStruct');
    chan_map_file = findChanMap(ks_dir);
    [unit_depth, ~, tip_depth, unit_x] = channelDepths(sessions{ii}, wv.max_site, data_dir, chan_map_file);
    keep = strcmp(wv.goodLabels(:), 'good') & ismember(wv.typeLabels(:), types) & ~isnan(unit_depth);
    depth = [depth; unit_depth(keep)]; %#ok<AGROW>
    shank = [shank; floor(unit_x(keep) / pair_pitch) + 1]; %#ok<AGROW>
    cell_type = [cell_type; wv.typeLabels(keep)]; %#ok<AGROW>
    fprintf('%s: tip %g um below the surface; good E %d, I %d, nonsomatic %d\n', sessions{ii}, tip_depth, ...
        nnz(keep & strcmp(wv.typeLabels(:), 'E')), nnz(keep & strcmp(wv.typeLabels(:), 'I')), nnz(keep & strcmp(wv.typeLabels(:), 'nonsomatic')))
end
m = load(chan_map_file);   % shank names: the front face's letter of each shank pair
map_pair = floor(m.xcoords(:) / pair_pitch) + 1;
n_shank = max(map_pair);
shank_name = arrayfun(@(p) char(unique(string(m.shank_letter(map_pair == p & string(m.side(:)) == "f")))), 1:n_shank, 'UniformOutput', false);
is_type = cellfun(@(t) strcmp(cell_type, t), types, 'UniformOutput', false);
session_str = strjoin(sessions, ', ');

%% 1. Each good unit at its shank and depth: units of one type at the same shank and depth (peak channel
% row) stack to the right in fixed steps, like a histogram; the three types are offset slightly in depth
x_step = 0.05;           % shank widths between stacked units
y_offset = [-4 0 4];     % um, E / I / nonsomatic
fig = figure('Name', 'Cell types vs depth'); ax = axes(fig); hold(ax, 'on')
h = gobjects(numel(types), 1);
for t = 1:numel(types)
    [~, ~, g] = unique([shank(is_type{t}), depth(is_type{t})], 'rows');   % units sharing shank, depth and type
    rank_in_stack = zeros(size(g)); for k = 1:numel(g); rank_in_stack(k) = nnz(g(1:k) == g(k)) - 1; end
    h(t) = scatter(ax, shank(is_type{t}) + rank_in_stack * x_step, depth(is_type{t}) + y_offset(t), 20, type_color{t}, 'LineWidth', 1);
end
set(ax, 'YDir', 'reverse', 'XTick', 1:n_shank, 'XTickLabel', shank_name, 'Box', 'off', 'TickDir', 'out', ...
    'FontSize', font_size, 'TitleFontSizeMultiplier', 1, 'LabelFontSizeMultiplier', 1)
xlim(ax, [0.8 n_shank + 0.8]); grid(ax, 'on')   % each stack starts at its shank's tick
xlabel(ax, 'shank (front face)'); ylabel(ax, 'depth below brain surface (\mum)')
legend(ax, h, compose('%s (%d)', string(types'), cellfun(@nnz, is_type')), 'Location', 'eastoutside', 'Box', 'off', 'FontSize', font_size)
title(ax, {['good units: ' session_str], ks_subdir}, 'Interpreter', 'none', 'FontWeight', 'normal')
drawnow

%% 2. Fraction of each type per depth bin: each shank, then pooled
edges = floor(min(depth) / depth_bin) * depth_bin : depth_bin : (floor(max(depth) / depth_bin) + 1) * depth_bin;
centers = edges(1:end-1) + depth_bin/2;
bin = discretize(depth, edges);
frac_offset = [-6 0 6];  % um, so the three types don't overlap
fig2 = figure('Name', 'Cell type fractions vs depth');
tl = tiledlayout(fig2, 1, n_shank + 1, 'TileSpacing', 'compact');
ax2 = gobjects(n_shank + 1, 1);
for p = 1:n_shank + 1
    if p <= n_shank; in = shank == p; name = shank_name{p}; else; in = true(size(shank)); name = 'all shanks'; end
    ax2(p) = nexttile(tl); hold(ax2(p), 'on')
    n_bin = accumarray(bin(in), 1, [numel(centers) 1]); has = n_bin > 0;
    for t = 1:numel(types)
        n_t = accumarray(bin(in), is_type{t}(in), [numel(centers) 1]);
        [frac, ci] = binofit(n_t(has), n_bin(has));
        errorbar(ax2(p), frac, centers(has)' + frac_offset(t), frac - ci(:, 1), ci(:, 2) - frac, 'horizontal', 'o', ...
            'Color', type_color{t}, 'CapSize', 3)
    end
    if p == n_shank + 1; text(ax2(p), ones(nnz(has), 1) * 1.02, centers(has), compose('n=%d', n_bin(has)), 'FontSize', font_size); end
    title(ax2(p), sprintf('%s (n = %d)', name, nnz(in)), 'FontWeight', 'normal')
    set(ax2(p), 'YDir', 'reverse', 'Box', 'off', 'TickDir', 'out', 'FontSize', font_size, ...
        'TitleFontSizeMultiplier', 1, 'LabelFontSizeMultiplier', 1)
    xlim(ax2(p), [0 1]); grid(ax2(p), 'on')
    if p > 1; ax2(p).YTickLabel = []; end
end
linkaxes(ax2, 'y'); ylim(ax2(1), edges([1 end]))
ylabel(ax2(1), 'depth below brain surface (\mum)'); xlabel(tl, 'fraction of good units (95% CI)', 'FontSize', font_size)
legend(ax2(end), compose('%s', string(types')), 'Location', 'southoutside', 'Box', 'off', 'FontSize', font_size)
title(tl, {['good units: ' session_str], sprintf('%s, %d um depth bins', ks_subdir, depth_bin)}, 'Interpreter', 'none', 'FontSize', font_size)
drawnow

%% Fractions per shank and pooled (all depths)
fprintf('\n%-12s %5s', 'shank', 'n'); fprintf('   %-24s', types{:}); fprintf('\n')
for p = 1:n_shank + 1
    if p <= n_shank; in = shank == p; name = shank_name{p}; else; in = true(size(shank)); name = 'all shanks'; end
    fprintf('%-12s %5d', name, nnz(in))
    for t = 1:numel(types)
        [frac, ci] = binofit(nnz(is_type{t}(in)), nnz(in));
        fprintf('   %3d  %.2f [%.2f-%.2f]   ', nnz(is_type{t}(in)), frac, ci(1), ci(2))
    end
    fprintf('\n')
end
