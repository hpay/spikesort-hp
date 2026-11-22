function wvStruct = getSessionWaveforms(dataDir, ksDir, only_good, raw, align, hp_pass)
% align (optional, default false): shift each unit's waveforms (all channels) so the
%   largest absolute extremum on the max channel, searched within the original window,
%   lands at sample round(spkOffset*fs)+1 (where most troughs already sit). Extra data
%   is extracted on both sides so the shifted window is filled with real data.
%   wvStruct.wf_shift = shift in samples (positive = extremum was late).
% raw (optional, for non-Intan data e.g. openEphys) struct with fields:
%   fbinary  - full path to the int16 raw binary
%   nCh      - total channels interleaved in the file (incl. ADC etc)
%   fs       - sample rate
%   chan_idx - 1-indexed ephys channels in the file to consider (e.g. 1:128).
%              waveFormsMean is restricted to these; max_site is still the
%              channel index in the raw file.
% If raw is omitted, reads Intan info.rhd + amplifier.dat in dataDir.
% hp_pass (optional): waveform high-pass cutoff (Hz), default 300 (see getSpikeWaveform)

% params
nWF = 1e3; % (HP modified from 1e3 to 1e2) % number of waveforms to get for each unit
spkOffset = 1.5e-3; % (SC modified from 1e-3 to 1.5e-3); % seconds before spike time to get
spkDur = 4e-3; % (s) total duration of waveform
dt = 0.5/1e3; % time step for ccg binning, SC modified from 1ms to 0.5ms

if exist('raw','var') && ~isempty(raw)
    nCh = raw.nCh;
    fs = raw.fs;
    fDat = raw.fbinary;
    chan_idx = raw.chan_idx(:)';
else
    % Get raw data info for Intan recording
    h = readIntanInfo(fullfile(dataDir,'info.rhd'));
    nCh = h.num_h.amplifier_channels;
    fs = h.sample_rate;
    fDat = fullfile(dataDir,'amplifier.dat');
    chan_idx = 1:nCh;
end

% set up files
spkSamp = readNPY(fullfile(ksDir, 'spike_times.npy'));
sID = readNPY(fullfile(ksDir,'spike_clusters.npy'));
[unit_ID,cluster_labels] = getPhyClusterLabels(ksDir);
if exist('only_good','var') && only_good
    ind = strcmp(cluster_labels,'good');
else
    ind = true(size(cluster_labels));
end
goodIDs = unit_ID(ind(:));
goodLabels = cluster_labels(ind(:));

% prepare variables
numUnits = length(goodIDs);
fprintf('getSessionWaveforms: extracting mean waveforms for %d units from %s\n', numUnits, fDat);
spkDurSamples = spkDur*fs; % total samples to pre-allocate for waveform
if ~exist('align','var') || isempty(align); align = false; end
if align
    n_pad = round((spkDur - spkOffset)*fs); % covers any shift of an extremum inside the original window
    i_target = round(spkOffset*fs) + 1;
else
    n_pad = 0;
end
wf_shift = zeros(numUnits,1);
waveFormsMean = nan(spkDurSamples,numel(chan_idx),numUnits);
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
    curUnitnSpikes = numel(curSpikeTimes); % numel: a 1-spike unit becomes a 1x0 row once its spike is removed
    if curUnitnSpikes == 0; continue; end % all spikes too close to the recording edges (waveform stays NaN)

    spikeTimesRP = curSpikeTimes(randperm(curUnitnSpikes));
    spikeTimeKeeps = double(sort(spikeTimesRP(1:min([nWF curUnitnSpikes]))));
    
    if ~exist('hp_pass', 'var'); hp_pass = []; end   % high-pass cutoff (Hz); [] = getSpikeWaveform's default, 300
    mean_wave_filt = getSpikeWaveform(mmf,spikeTimeKeeps, fs, spkOffset + n_pad/fs, spkDur + 2*n_pad/fs, hp_pass);
    mean_wave_filt = mean_wave_filt(:,chan_idx);
    orig_win = n_pad + (1:spkDurSamples);

    % Get channel with largest amplitude (in the original window), take that as the waveform
    amps = max(mean_wave_filt(orig_win,:))-min(mean_wave_filt(orig_win,:));
    [max_val,i_max] = max(amps);
    if align
        [~, i_ext] = max(abs(mean_wave_filt(orig_win,i_max)));
        wf_shift(thisUnit) = i_ext - i_target;
        orig_win = orig_win + wf_shift(thisUnit);
    end
    waveFormsMean(:,:,thisUnit) = mean_wave_filt(orig_win,:);
    max_site(thisUnit) = chan_idx(i_max); % Max site is the raw file channel number, so 1 = Intan A-000 / openEphys CH1
    mxWF(thisUnit,:) = waveFormsMean(:,i_max,thisUnit)';
    
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
if align; wvStruct.wf_shift = wf_shift; end
end