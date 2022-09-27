function results = get_session_waveforms(dataDir, ksDir, only_good)
% params
nWF = 1e3; % (HP modified from 1e3 to 1e2) % number of waveforms to get for each unit
spkOffset = 1.5e-3; % (SC modified from 1e-3 to 1.5e-3); % seconds before spike time to get
spkDur = 4e-3; % seconds after spike time to get
dt = 0.5/1e3; % time step for ccg binning, SC modified from 1ms to 0.5ms

% Get raw data info for Intan recording
h = readIntanInfo(fullfile(dataDir,'info.rhd'));
nCh = h.num_h.amplifier_channels;
fs = h.sample_rate;

% set up files
fDat = fullfile(dataDir,'amplifier.dat');
spkSamp = readNPY(fullfile(ksDir, 'spike_times.npy'));
sID = readNPY(fullfile(ksDir,'spike_clusters.npy'));
[cIDs,cluster_labels] = get_phy_cluster_labels(ksDir);
if only_good
    ind = cluster_labels==2;
else
    %         ind = cluster_labels>0 & cluster_labels<3;
    ind = isfinite(cluster_labels);
end
goodIDs = cIDs(ind(:));
goodLabels = cluster_labels(ind(:));

% prepare variables
numUnits = length(goodIDs);
spkDurSamples = spkDur*fs; % total samples to pre-allocate for waveform
waveFormsMean = nan(spkDurSamples,nCh,numUnits);
meanRate = nan(numUnits,1);
medISI = nan(numUnits,1);
contam = nan(numUnits,1); % Why are we calculating contam here? HP
filenamestruct = dir(fDat);
nSamp = filenamestruct.bytes/(nCh*2);  % Number of samples per channel (int16 is 2 bytes each sample)
mmf = memmapfile(fDat, 'Format', {'int16', [nCh nSamp], 'x'});

maxTime = double(max(spkSamp))/fs;

% read in spikes for all units
for thisUnit=1:numUnits
    curUnitID = goodIDs(thisUnit);
    curSpikeTimes = double(spkSamp(sID==curUnitID))/fs;
    meanRate(thisUnit) = length(curSpikeTimes)./maxTime;
    medISI(thisUnit) = median(diff(curSpikeTimes));
    if ~isempty(curSpikeTimes)
        [K, Qi, Q00, Q01, rir] = ccg(curSpikeTimes, curSpikeTimes, 500, dt);
        contam(thisUnit) = min(Qi/(max(Q00, Q01)));
        
        % Exclude spikes whos waveform extends beyond recording
        curSpikeTimes((curSpikeTimes+spkDur*10)*fs > nSamp) = []; % spikes extending past recording
        curSpikeTimes((curSpikeTimes-spkOffset*10)*fs < 1) = []; % spikes starting before recording
        curUnitnSpikes = size(curSpikeTimes,1);
        
        spikeTimesRP = curSpikeTimes(randperm(curUnitnSpikes));
        spikeTimeKeeps = double(sort(spikeTimesRP(1:min([nWF curUnitnSpikes]))));
        
        mean_wave_filt = getSpikeWaveform(mmf, ...
            spikeTimeKeeps, fs, spkOffset, spkDur);
        waveFormsMean(:,:,thisUnit) = mean_wave_filt;
    end
    disp(['Completed ' int2str(thisUnit) ' units of ' int2str(numUnits) '.']);
end

% Get top PC of waveforms across channels
pcWF = nan(numUnits,spkDurSamples);

for i=1:numUnits
    [e,s,l] = pca(waveFormsMean(:,:,i)');
    % correct sign so that largest value is negative
    if max(e(:,1)) > min(e(:,1))
        wvSign = -1;
    else
        wvSign = 1;
    end
    pcWF(i,:) = wvSign * e(:,1)';
end

% Get channel with largest amplitude, take that as the waveform
amps = squeeze(max(waveFormsMean)-min(waveFormsMean));
[max_val,max_site] = max(amps);
max_site = max_site(:); % Max site is the Intan channel number, so 1 = A-000
mxWF = nan(numUnits,spkDurSamples);
for i=1:numUnits
    mxWF(i,:) = waveFormsMean(:,max_site(i),i)';
end


results = struct('mxWF',mxWF,'max_site',max_site,'pcWF',pcWF,'meanRate',meanRate,...
    'waveFormsMean',waveFormsMean,'spkDur',spkDur,'spkOffset',spkOffset,'goodIDs',goodIDs,...
    'goodLabels',goodLabels,'medISI',medISI,'contam',contam);
end