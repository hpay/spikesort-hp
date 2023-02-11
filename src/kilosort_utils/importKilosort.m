function S = importKilosort(ksDir, fs, option_only_good, option_orig, session_label)
% IMPORTKILOSORT import sorted spike times and waveforms and convert to dat
% structure
%
% S = IMPORTKILOSORT(filepath, fs) imports all available channels of sorted data
% 	from the folder filepath, recorded at samplerate fs (usually 30000)
% option_only_good: Set 1 to only load "good" cells (post manual phy sorting, if available)
% option_orig: ignore any phy sorting and use categories assigned by kilosort
%
% session_label: a prefix added to each cell name (before K000 assigned unit number etc)
%
% Example:
% S = importKilosort('Z:\Hannah\ephys\project2\HC05_220825\kilosort2_output', 3e4, 1)

% read in spike times and get IDs and kilosort labels
tm = readNPY(fullfile(ksDir,'spike_times.npy'));
tm = double(tm)/fs;
sID = readNPY(fullfile(ksDir,'spike_clusters.npy'));
if exist('option_orig','var') && option_orig
    [cIDs, cluster_labels] = get_kilo_cluster_labels(ksDir); % 0 = noise, 1 = 'mua', 2 = 'good'
else
    [cIDs, cluster_labels] = get_phy_cluster_labels(ksDir); % 0 = noise, 1 = 'mua', 2 = 'good'
end

% Only look at "good" clusters -- if sorted with Phy, this will be the
% final sorting
if option_only_good
    cIDs = cIDs(strcmp(cluster_labels,'good'));
    cluster_labels = cluster_labels(strcmp(cluster_labels,'good'));
else % Discard "noise" sorted clusters
    cIDs = cIDs(~strcmp(cluster_labels,'noise'));
    cluster_labels = cluster_labels(~strcmp(cluster_labels,'noise'));
end

% Load these from waveformStruct.mat!
waveform_file = fullfile(ksDir, 'waveformStruct.mat');
if exist(waveform_file,'file')
    wvStruct = getfield(load(waveform_file),'wvStruct');
else
    wvStruct = [];
    wvStruct.max_site = NaN(length(cIDs),1);
    wvStruct.mxWF = NaN(length(cIDs),1);
    wvStruct.spkOffset = NaN;
end

% Loop through units.
S = dat;
for ii = 1:length(cIDs)
    
    % Cell ID
    cID = cIDs(ii);
    chan_val = wvStruct.max_site(ii); % intan channel index of largest waveform. So 1 = A-000.
    
    % Get spike times from this unit
    tt = tm(sID==cID);
    
    % Store spikes
    if exist('session_label','var')&& ~isempty(session_label)
        chan_label =  sprintf('%s_K%03i',session_label, cID); % To indicate KS assigned unit number
    else
        chan_label =  sprintf('K%03i', cID); % To indicate KS assigned unit number
    end
    tstart = min(tt);
    tend = max(tt);
    
    % Store as dat structure
    %  S(ii) = dat(tt, chan_label, chan_val, 'event', tstart, tend, 's') ;
    
    % Include waveform and store as dat structure
    %  waveform = wvStruct.pcWF(ii,:); % pcWF stores top PC (NOT in uV)!
    %  info.waveunit = 'A.U.';
    waveform = wvStruct.mxWF(ii,:); % mxWF stores largest waveform!
    info.waveunit = 'uV';
    info.wavefreq = fs;
    info.prethresh = fs*wvStruct.spkOffset;
    S(ii) = dat(tt, chan_label, chan_val, 'event', tstart, tend, 's', cID, waveform, info) ;
    
end


