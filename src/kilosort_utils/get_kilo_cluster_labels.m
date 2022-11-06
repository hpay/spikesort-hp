function [cIDs,cluster_labels] = get_kilo_cluster_labels(ksDir)
% equivalent to reading cluster_group.tsv from kilosort to get ORIGINAL cluster labels
% according to the following code:
% 1 - 'mua'
% 2 - 'good'/well-isolated

fid = fopen(fullfile(ksDir,'cluster_KSlabel.tsv'));
C = textscan(fid, '%s%s');
fclose(fid);

cIDs = cellfun(@str2num, C{1}(2:end), 'uni', false);
ise = cellfun(@isempty, cIDs);
cIDs = [cIDs{~ise}];

isMUA = cellfun(@(x)strcmp(x,'mua'),C{2}(2:end));
isGood = cellfun(@(x)strcmp(x,'good'),C{2}(2:end));
cluster_labels = zeros(size(cIDs));

cluster_labels(isMUA) = 1;
cluster_labels(isGood) = 2;
