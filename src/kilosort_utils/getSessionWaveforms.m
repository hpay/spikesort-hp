function wvStruct = getSessionWaveforms(dataDir, ksDir, only_good)
% params
nWF = 1e3; % (HP modified from 1e3 to 1e2) % number of waveforms to get for each unit
spkOffset = 1.5e-3; % (SC modified from 1e-3 to 1.5e-3); % seconds before spike time to get
spkDur = 4e-3; % (s) total duration of waveform
dt = 0.5/1e3; % time step for ccg binning, SC modified from 1ms to 0.5ms

% Get raw data info for Intan recording
h = readIntanInfo(fullfile(dataDir,'info.rhd'));
nCh = h.num_h.amplifier_channels;
fs = h.sample_rate;

% set up files
fDat = fullfile(dataDir,'amplifier.dat');
spkSamp = readNPY(fullfile(ksDir, 'spike_times.npy'));
sID = readNPY(fullfile(ksDir,'spike_clusters.npy'));
[unit_ID,cluster_labels] = getPhyClusterLabels(ksDir);
if only_good
    ind = strcmp(cluster_labels,'good');
else
    ind = true(size(cluster_labels));
end
goodIDs = unit_ID(ind(:));
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
max_site = NaN(numUnits,1);
mxWF = nan(numUnits,spkDurSamples);
nSpikes = NaN(numUnits,1);

% read in spikes for all units
for thisUnit=1:numUnits
    curUnitID = goodIDs(thisUnit);
    curSpikeTimes = double(spkSamp(sID==curUnitID))/fs;
    meanRate(thisUnit) = length(curSpikeTimes)./maxTime;
    medISI(thisUnit) = median(diff(curSpikeTimes));
    nSpikes(thisUnit) = length(curSpikeTimes);
    if isempty(curSpikeTimes); continue; end
    
    [K, Qi, Q00, Q01, Ri] = ccg(curSpikeTimes, curSpikeTimes, 500, dt);
    contam(thisUnit) = min(Qi/(max(Q00, Q01)));
    % Q00:  measure of outer shoulder height, norm Poisson expectation
    % Q01: inner shoulder height, norm. by Poisson expectation
    % Qi: measure of inner refractory period height, different window
    % sizes (1,2,3 bins etc)
        
    % Exclude spikes whos waveform extends beyond recording
    curSpikeTimes((curSpikeTimes+spkDur*10)*fs > nSamp) = []; % spikes extending past recording
    curSpikeTimes((curSpikeTimes-spkOffset*10)*fs < 1) = []; % spikes starting before recording
    curUnitnSpikes = size(curSpikeTimes,1);
    
    spikeTimesRP = curSpikeTimes(randperm(curUnitnSpikes));
    spikeTimeKeeps = double(sort(spikeTimesRP(1:min([nWF curUnitnSpikes]))));
    
    mean_wave_filt = getSpikeWaveform(mmf,spikeTimeKeeps, fs, spkOffset, spkDur);
    waveFormsMean(:,:,thisUnit) = mean_wave_filt;
    
    % Get channel with largest amplitude, take that as the waveform
    amps = max(mean_wave_filt)-min(mean_wave_filt);
    [max_val,max_site(thisUnit)] = max(amps); % Max site is the Intan channel number, so 1 = A-000
    mxWF(thisUnit,:) = waveFormsMean(:,max_site(thisUnit),thisUnit)';
    
    %         K((length(K)+1)/2) = 0;
    %         figure; subplot(2,2,1); plot(K); subplot(2,2,2); plot(Qi); title('Qi'); subplot(2,2,3); plot(Ri); title('Ri')
    %         subplot(2,2,4); plot(mxWF(thisUnit,:),'k')
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

wvStruct = struct('mxWF',mxWF,'max_site',max_site,'pcWF',pcWF,'meanRate',meanRate,...
    'waveFormsMean',waveFormsMean,'spkDur',spkDur,'spkOffset',spkOffset,'goodIDs',goodIDs,...
    'goodLabels',{goodLabels},'medISI',medISI,'contam',contam,'nSpikes',nSpikes);
end