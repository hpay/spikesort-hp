function [stats, goodIDs, mxWF, tWF] = getGMM_stats(all_dir, option_only_good)


if istable(all_dir)
    R = all_dir;
    meanRate = R.rate;
    mxWF = cell2mat([R.waveform]')';
    goodIDs = ones(size(meanRate));
    fs = 1e6;
    tWF = ((0:size(mxWF,2)-1)/fs - 3e-4)*1e3;
else
    
    if ~iscell(all_dir); all_dir = {all_dir}; end
    for ii = 1:length(all_dir)
        ks_dir = fullfile(all_dir{ii},'kilosort2_output');
        wvStruct = getfield(load(fullfile(ks_dir,'waveformStruct.mat')),'wvStruct');
        wvStruct.goodIDs = wvStruct.goodIDs(:);
        wvStruct.goodLabels = wvStruct.goodLabels(:);
        
        % Only keep "good" clusters. If sorted with Phy, this will be the final sorting
        if option_only_good
            [unit_ID, cluster_labels] = get_phy_cluster_labels(ks_dir); % 0 = noise, 1 = 'mua', 2 = 'good'
            if length(unit_ID) ~= length(wvStruct.goodIDs); error('Run waveform extraction after manual curation'); end
            mask = strcmp(cluster_labels,'good');
            wvStruct.mxWF = wvStruct.mxWF(mask,:);
            wvStruct.max_site = wvStruct.max_site(mask,:);
            wvStruct.pcWF = wvStruct.pcWF(mask,:);
            wvStruct.meanRate = wvStruct.meanRate(mask,:);
            wvStruct.waveFormsMean = wvStruct.waveFormsMean(:,:,mask);
            wvStruct.goodIDs = wvStruct.goodIDs(mask);
            wvStruct.goodLabels = wvStruct.goodLabels(mask);
            wvStruct.medISI = wvStruct.medISI(mask);
            wvStruct.contam = wvStruct.contam(mask);
        end
        
        % Add mxWF_d
        theseWaveforms(ii) = wvStruct;
    end
    fs = size(wvStruct.mxWF,2)/wvStruct.spkDur;
    
    tWF = linspace(-theseWaveforms(1).spkOffset,-theseWaveforms(1).spkOffset+theseWaveforms(1).spkDur,...
        fs*theseWaveforms(1).spkDur);
    
    meanRate = cat(1,theseWaveforms.meanRate);
    mxWF = cat(1,theseWaveforms.mxWF);
    goodIDs = cat(1,theseWaveforms.goodIDs);
    
end

%% Compute stats
mxWF_d = diff(mxWF')';
numUnits = length(meanRate);
stats = array2table(NaN(numUnits,4)); stats.Properties.VariableNames = {'rate_log10','duration','asymmetry','pk_deriv'};
for nUnit = 1:numUnits
    % get trough
    [tr, trInd] = min(mxWF(nUnit,:));
    
    % get peak (after trough)
    [pk, pkInd] = max(mxWF(nUnit, trInd:end));
    pkInd = pkInd + trInd - 1;
    
    % get pre-peak (before trough)
    pk_pre = max(mxWF(nUnit, 1:trInd));
    
    % get derivative trough and peak
    [tr_d,trInd_d] = min(mxWF_d(nUnit,:));
    [pk_d,pkInd_d] = max(mxWF_d(nUnit,:));
    
    % concatenate the following statS
    stats.rate_log10(nUnit) = log10(meanRate(nUnit));
    stats.duration(nUnit) =  pkInd-trInd;               % duration
    stats.asymmetry(nUnit) = (pk-pk_pre)./(pk+pk_pre);  % asymmetry
    stats.pk_deriv(nUnit) = log10(abs(pk_d./tr_d));     % pk derivative ratio
end


end
