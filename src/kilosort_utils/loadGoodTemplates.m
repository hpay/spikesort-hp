function [good, wf, pos, peak_ch, amp_max, tmpl, group] = loadGoodTemplates(ks_dir, groups)
% LOADGOODTEMPLATES  Unwhitened templates of the good units (cluster_group.tsv) of a Kilosort 2 or 4 sort.
%   groups   (optional, default {'good'}) cluster_group.tsv labels to include, e.g. {'good', 'mua'}
%   good     cluster ids of the units included; group = their labels
%   wf       units x samples x channels, the template most of each unit's spikes came from, times
%            whitening_mat_inv.npy
%   pos      channels x [x y] (channel_positions.npy)
%   peak_ch  channel with the largest peak-to-peak template, per unit; amp_max = that peak-to-peak
%   tmpl     0-based template index of each unit (rows of similar_templates.npy)

templates = readNPY(fullfile(ks_dir, 'templates.npy'));            % templates x samples x channels (whitened)
w_inv = readNPY(fullfile(ks_dir, 'whitening_mat_inv.npy'));
pos = readNPY(fullfile(ks_dir, 'channel_positions.npy'));
spike_templates = readNPY(fullfile(ks_dir, 'spike_templates.npy'));
spike_clusters = readNPY(fullfile(ks_dir, 'spike_clusters.npy'));
T = readtable(fullfile(ks_dir, 'cluster_group.tsv'), 'FileType', 'text', 'Delimiter', '\t');
if ~exist('groups', 'var'); groups = {'good'}; end
keep = ismember(T.group, groups);
good = T.cluster_id(keep); group = T.group(keep);

[n_samp, n_ch] = deal(size(templates, 2), size(templates, 3));
wf = zeros(numel(good), n_samp, n_ch);
tmpl = zeros(numel(good), 1);
for k = 1:numel(good)
    tmpl(k) = mode(spike_templates(spike_clusters == good(k)));
    wf(k, :, :) = squeeze(templates(tmpl(k) + 1, :, :)) * w_inv;
end
amp = squeeze(max(wf, [], 2) - min(wf, [], 2));                    % units x channels, peak to peak
[amp_max, peak_ch] = max(amp, [], 2);
end
