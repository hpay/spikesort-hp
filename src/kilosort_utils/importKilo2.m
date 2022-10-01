function S = importKilo2(ksDir, fs, option_only_good, session_label)
% IMPORTKILO2 import sorted spike times and waveforms and convert to dat
% structure
%
% S = importKilo2(filepath) imports all available channels of sorted data
% from the folder filepath.
%
% Example:
% S = importKilo2('Z:\Hannah\ephys\project2\HC05_220825\raw_220825_131723\kilo2.0', 3e4, 1)

% read in spike times and get IDs and kilosort labels
tm = readNPY(fullfile(ksDir,'spike_times.npy'));
tm = double(tm)/fs;
sID = readNPY(fullfile(ksDir,'spike_clusters.npy'));
[cIDs, cluster_labels] = get_phy_cluster_labels(ksDir); % 0 = noise, 1 = 'mua', 2 = 'good'

% Only look at "good" clusters?
if option_only_good
    good_label = 2;
    cIDs = cIDs(cluster_labels==good_label);
    cluster_labels = cluster_labels(cluster_labels==good_label);
end

% Get waveforms - SLOW! Now doing this separately right after spike sorting
% wvStruct = getSessionWaveforms(dataDir, ksDir, option_only_good);

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

S = dat;
% Loop through units.
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


