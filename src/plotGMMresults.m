function plotGMMresults(tWF, mxWF, stats, id_text, good_all, session_all, Ranatomy, channel_all, unitID_all)
%% Plot results


% TODO: change and id where "good" label = noise to unknown
[id_label, ~, id_num] = unique(id_text);
n = height(stats);
fs = 1/mean(diff(tWF));

cE = [182 47 100]/255;  % pink
cI = [29,29,154]/255;   % blue
cNS = .6*[1 1 1];


disp('number "E" or "I" cells classified as "noise":')
disp(nnz(ismember(id_text,{'E','I'}) & ismember(good_all,'noise')))

ds = 1; % Downsample all cells so plot is interpretable
mask_ds = false(n,1); mask_ds(1:ds:n) = true;
mask_singleunit = ~ismember(id_text,{'extreme_outlier','unknown'}) & strcmp(good_all,'good');
mask_E = strcmp(id_text,'E');
mask_I = strcmp(id_text,'I');
mask_intermediate = strcmp(id_text,'intermediate');
mask_outlier = strcmp(id_text,'moderate_outlier');
mxWFnorm = mxWF./-min(mxWF,[],2);


figure;
msize  = 5;
malpha = .5;
mask = mask_singleunit & ~(mask_E|mask_I) & mask_ds;
scatter3(stats.duration(mask)/fs*1e3, stats.asymmetry(mask),10.^stats.rate_log10(mask),msize,cNS,'MarkerEdgeAlpha',malpha); hold on
mask = mask_singleunit & mask_E & mask_ds;
scatter3(stats.duration(mask)/fs*1e3, stats.asymmetry(mask),10.^stats.rate_log10(mask),msize,cE,'MarkerEdgeAlpha',malpha); hold on
mask = mask_singleunit & mask_I & mask_ds;
scatter3(stats.duration(mask)/fs*1e3, stats.asymmetry(mask),10.^stats.rate_log10(mask),msize,cI,'MarkerEdgeAlpha',malpha); hold on
xlabel('Duration tr-pk(ms)'); ylabel('Asymmetry'); zlabel('Firing rate (Hz)')
set(gca,'ZScale','log')
xlim([0  1]); ylim([-1 1]); zlim(10.^[-2.2 1.8]); 
set(gca,'XTick',0:.5:1,'YTick',-1:.5:1,'ZTick',10.^[-2:1])
set(gca,'ZMinorGrid','off');
axis vis3d
set(gca,'YDir','reverse')
az = 45;
el = 20;
view(az, el)


%% TEMP select interesting points in 2D
%{
R_z = [cosd(az), -sind(az), 0;
    sind(az), cosd(az), 0;
    0, 0, 1];
R_x = [1, 0, 0;
    0, cosd(el), -sind(el);
    0, sind(el), cosd(el)];
R = R_z * R_x;
mask = mask_singleunit & (mask_E|mask_I); %& mask_ds;
points3d = [stats.duration/fs*1e3, stats.asymmetry, stats.rate_log10];
points3d(~mask,:) = NaN;

rotatedPoints3D = (R * points3d')';
xx = rotatedPoints3D(:,1);
yy = rotatedPoints3D(:,3);

figure;
subplot(1,2,1)
scatter(xx,yy,msize,'k'); hold on;
title('Draw Polygon to Select Points');

% Draw polygon
for ii = 1
    subplot(1,2,1)
    h = impoly; position = wait(h); % Wait until the polygon is closed
    colorOrder = get(gca, 'ColorOrder');
    cc = colorOrder(mod(get(gca, 'ColorOrderIndex')-1, size(colorOrder, 1)) + 1, :);
    
    % Find points inside the polygon
    mask_select = inpolygon(xx, yy, position(:,1), position(:,2)) & mask;
    idx_select = find(mask_select);
    scatter(xx(mask_select), yy(mask_select),msize, cc)
    subplot(122)
    plot(tWF, mxWFnorm(idx_select(1:10:end),:),'Color',cc,'LineWidth',.5); hold on
    fixticks; drawnow
end
%}


%% Plot session name v. cell depth highlighting selected cells
% figure
% [session_labels, ~, session_id] = unique(session_all);
% mask = mask_singleunit & (mask_E|mask_I) & ~mask_select;
% scatter(session_id(mask), Ranatomy.DV(mask),20,'k'); hold on
% mask = mask_singleunit & (mask_E|mask_I) & mask_select;
% scatter(session_id(mask)+.1, Ranatomy.DV(mask),20,'r'); hold on
% set(gca,'XTick',1:length(session_labels), 'XTickLabel',session_labels,'FontSize', 7, 'TickLabelInterpreter', 'none');
% ylabel('DV (mm)')
% xtickangle(90);
% legend('other cells','strange E region')


%% Print out bad cells from example session
% mask = find(mask_select& strcmp(session_all, 'HC18R_231213'));
% unitID_all(mask) % starts at 0
% channel_all(mask) % starts at 1
% Ranatomy.DV(mask)
% good_all(mask)


%% Plot waveforms
%{
rng(1)
nplot = 30;

figure
ss = 1;
hs = subplot(2,2,ss);
mask = find(mask_singleunit & mask_E);
mask = mask(randperm(numel(mask), nplot),:);
cc = cE;
 plot(tWF,mxWFnorm(mask, :)','Color',cc);   hold on,
 plot(tWF,mean(mxWFnorm(mask, :)),'Color','k','LineWidth',2);
axis square; title('Excitatory'); ylabel('norm. amplitude');


ss = ss+1;
hs(ss) = subplot(2,2,ss);
mask = find(mask_singleunit & mask_I);
mask = mask(randperm(numel(mask), nplot),:);
cc = cI;
 plot(tWF,mxWFnorm(mask, :)','Color',cc);   hold on,
 plot(tWF,mean(mxWFnorm(mask, :)),'Color','k','LineWidth',2);
axis square; title('Inhibitory'); ylabel('norm. amplitude');

ss = ss+1;
hs(ss) = subplot(2,2,ss);

plot(tWF,zscore(mxWF(id_num==0, :)'),'Color','k');  hold on; 
plot(tWF,zscore(mxWF(id_num==3, :)'),'Color',[ 0.6000    0.4000    1.0000]); hold on;
axis square; title('Intermediate/unknown')

ss = ss+1;
hs(ss) = subplot(2,2,ss);
 plot(tWF,zscore(mxWF(id_num==4, :)'),'Color',[0.8 .2 1]); hold on; 
     plot(tWF,zscore(mxWF(id_num==5, :)'),'Color','k');  
%     legend('intermediate', 'moderate_outlier','extreme_outlier','unknown')
axis square; title('moderate outlier/extreme outlier')
linkaxes(hs); xlim(tWF([1 end]));
ylim([-6 3])
grid(hs,'on')
drawnow

%}