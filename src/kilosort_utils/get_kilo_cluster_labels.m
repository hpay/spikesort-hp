function [unit_ID,cluster_labels] = get_kilo_cluster_labels(ksDir)
% equivalent to reading cluster_group.tsv from kilosort to get ORIGINAL cluster labels

fid = fopen(fullfile(ksDir,'cluster_KSlabel.tsv'));
C = textscan(fid, '%s%s');
fclose(fid);

unit_ID = cellfun(@str2num, C{1}(2:end), 'uni', false);
ise = cellfun(@isempty, unit_ID);
unit_ID = [unit_ID{~ise}];
cluster_labels = C{2}(2:end);