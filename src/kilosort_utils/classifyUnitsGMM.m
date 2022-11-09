function [idx,nlogL,P,logpdf] = classifyUnitsGMM(wvStruct)

mxWF = cat(1,wvStruct.mxWF);
mxWF_d = cat(1,wvStruct.mxWF_d);
meanRate = cat(1,wvStruct.meanRate);
medISI = cat(1,wvStruct.medISI);
unitLabel = cat(2,wvStruct.goodLabels)';

% get stats
numUnits = length(meanRate);
allStats = nan(numUnits, 5);
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
    
    % concatenate the following stats
    %(1) pk derivative ratio (2) duration (3) derivative duration (4)
    %asymmetry (5) Burstiness (inverse-median ISI minus meanRate)
    allStats(nUnit,:) = [log10(abs(pk_d./tr_d)),pkInd-trInd,pkInd_d-trInd_d,...
    (pk-pk_pre)./(pk+pk_pre),log10(1./medISI(nUnit)) - log10(meanRate(nUnit))];
end

load('latestGMM_cellType.mat','gm'),
[idx,nlogL,P,logpdf] = gm.cluster(allStats);
