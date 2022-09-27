% Copied from Selmaan, haven't done anything really with this yet

close all; clear;

% directories
baseDir = dir('Z:\Hannah\ephys\project2\HC05_220816\raw*');
baseDir = fullfile(baseDir.folder,baseDir.name);
% ksDir = fullfile(baseDir,'kilo2.0'); % Selmaan's organization - to implement
ksDir = baseDir;

% read in spike times and get IDs and kilosort labels
tm = readNPY(fullfile(ksDir,'spike_times.npy'));
tm = double(tm)/30e3;
sID = readNPY(fullfile(ksDir,'spike_clusters.npy'));
[cIDs,cluster_labels] = get_phy_cluster_labels(ksDir); % 1 - 'mua'; 2 - 'good'

% classify waveforms using GMM
wvStruct = get_session_waveforms(baseDir, ksDir, 0);
[idx,~,P,logpdf] = classifyUnitsGMM(wvStruct);

% idx(logpdf<-10) = 4; % heuristic to create a new group for all highly unusual waveforms

% get rid of (near-)empty clusters
emptyNeur = find(wvStruct.meanRate < 0.01);
cIDs(emptyNeur) = [];
cluster_labels(emptyNeur) = [];
idx(emptyNeur) = [];
indI = find(idx==2);
indE = find(idx==1);
contam = wvStruct.contam(setdiff(1:length(wvStruct.contam),emptyNeur));
%%
f = fopen(fullfile(baseDir,'digitalin.dat'));
frameTrig = fread(f,'*uint16');
fclose(f);
frameTrig = frameTrig > 1/2;
frameOnsets = find(frameTrig(2:end) & ~frameTrig(1:end-1));
frameOnsets = frameOnsets/30e3;
invalidIFI = diff(frameOnsets)>1.1/60;
if any(invalidIFI)
    error('Check IFIs!'),
end

%%
xBeh = [frameOnsets; frameOnsets(end)+1/60];
spks = nan(length(xBeh)-1,length(cIDs));
for i=1:length(cIDs)
    theseSpks = tm(sID==cIDs(i));
    spks(:,i) = histcounts(theseSpks, xBeh);
end


%% count active neurons in 50ms bins
goodNeur = find(contam<1/3);
goodI = intersect(goodNeur,indI);
goodE = intersect(goodNeur,indE);

actNeur = sum(spks(:,goodNeur)>0,2);
actNeur = movmean(actNeur,3)*3;
tmp = spks(:,goodNeur)>0;
for i=1:size(tmp,2)
    tmp(:,i) = circshift(tmp(:,i),randi(size(tmp,1)));
end
actShuf = sum(tmp,2);
actShuf = movmean(actShuf,3)*3;

figure,plot(actNeur),
shufBounds = prctile(actShuf,[1,99]);
hold on,plot(xlim,[1,1]*shufBounds(1),'r--','LineWidth',2),
plot(xlim,[1,1]*shufBounds(2),'r--','LineWidth',2),
figure,histogram(actNeur,'Normalization','pdf','BinMethod',"integers")
hold on,histogram(actShuf,'Normalization','pdf','BinMethod',"integers")
set(gca,'YScale',"log")

figure,hold on,
plot(1e3*(-60:60)/60,xcov(sum(spks(:,goodI)>0,2),60,'coef'))
plot(1e3*(-60:60)/60,xcov(sum(spks(:,goodE)>0,2),60,'coef'))
plot(1e3*(-60:60)/60,xcov(sum(spks(:,goodE)>0,2),sum(spks(:,goodI)>0,2),60,'coef'))
plot(1e3*(-60:60)/60,xcov(sum(spks(:,goodNeur)>0,2),60,'coef'))


%%
goodNeur = find(contam<1/3);
rCoh = nan(length(cIDs),1);
% normSpks = zscore(spks);
normSpks = spks>0;
for i=1:length(goodNeur)
    thisNeur = goodNeur(i);
    otherNeur = setdiff(goodNeur,thisNeur);
    thisCount = mean(normSpks(:,otherNeur),2);
    rCoh(goodNeur(i)) = corr(thisCount, spks(:,thisNeur));
end

figure,semilogx(mean(spks(:,goodNeur))*60,rCoh(goodNeur),'.')

%%
% goodNeur = find(alignedData.contam<1/3);
goodNeur = intersect(find(alignedData.contam<1/2),alignedData.indE);
hd = mean(alignedData.smPts(:,1:2,[9,14]),3);
iMax = size(hd,1);
f = 15*25.4;
cd('C:\Users\Selmaan\Documents\MATLAB\scratch code')
for i=1:length(goodNeur)
    iSpks = repelem(1:iMax, alignedData.spks(1:iMax, goodNeur(i)))';
    [map, info, bin_centers, pos_time_smooth] = placeCellAnalysis(...
        hd(:,1)*f, hd(:,2)*f, 60, iSpks, 0, [-.9, .9]*f, 40, 13, 0.25, 1);
    title(sprintf('Cell #%d, Rate %0.2f, Info %0.2f',...
        goodNeur(i),length(iSpks)/(iMax/60), info)),
    pause,
end

%%

gH = bk(:,3)-hd(:,3);
gD = atan2d(bk(:,2)-hd(:,2),bk(:,1)-hd(:,1));
nBins = 100;
xBins = linspace(-180,180,nBins);
binWidth = 10;
goodNeur = find(contam<1/3);
normSpks = zscore(spks(:,goodNeur)); %spks(:,goodNeur)./mean(spks(:,goodNeur));
tmp = nan(nBins, length(goodNeur));
for i=1:nBins
    ind = (gD > xBins(i)-binWidth) & (gD < xBins(i)+binWidth);
    tmp(i,:) = mean(normSpks(ind,:));
end
figure,plot(xBins,tmp)

%%
nCoef = 7;
posTarg = hd(:,1:2);
posCenters = linspace(-.9,.9,nCoef);
posWidth = median(diff(posCenters));
[xPos, yPos] = meshgrid(posCenters, posCenters);
X = nan(size(posTarg,1), nCoef^2);
for x=1:nCoef
    for y=1:nCoef
        thisCent = [xPos(x,y), yPos(x,y)];
        thisFunc = @(pos) exp(-sum((pos-thisCent).^2,2) ./ posWidth.^2);
        X(:,x + (y-1)*nCoef) = thisFunc(posTarg);
    end
end

%%
smooth_xx = movmean(hd(:,1),61);
smooth_yy = movmean(hd(:,2),61);
speed = [sqrt(diff(smooth_xx).^2+diff(smooth_yy).^2)*60; nan]*15*25.4;
speedInd = movmean(speed,301) > 25;

normSpks = zscore(spks(1:iMax,goodNeur));

%%
nNeur = 106;
nDecimate = 10;
nPts = size(X,1)-mod(size(X,1),nDecimate);
inclMask = sum(X) > 1e2;

thisX = squeeze(mean(reshape(X(1:nPts,:), nDecimate, nPts/nDecimate, size(X,2))));
thisX = thisX(:,inclMask);
thisY = squeeze(sum(reshape(spks(1:nPts,nNeur), nDecimate, nPts/nDecimate)))';

mdl = fitglm(thisX, thisY, 'Distribution',  'poisson');
y_hat = mdl.predict(X(:,inclMask));
figure,scatter(posTarg(:,1),posTarg(:,2), 10, y_hat/nDecimate, 'filled'),