function [unit_ID,cluster_labels] = getPhyClusterLabels(ksDir)
% reads cluster_group.tsv from kilosort/phy, and gets cluster labels

fid = fopen(fullfile(ksDir,'cluster_group.tsv'));
C = textscan(fid, '%s%s');
fclose(fid);

unit_ID = cellfun(@str2num, C{1}(2:end), 'uni', false);
ise = cellfun(@isempty, unit_ID);
unit_ID = [unit_ID{~ise}];
unit_ID = unit_ID(:);
cluster_labels = C{2}(2:end);