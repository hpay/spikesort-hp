function plotBasicResults(T, data_dir, ks_dir_fun)
% T.date
% T.surgery
% T.bird
% T.depth
% T.manually_sorted




%% Colors
birds = unique(T.bird);
cs = turbo(length(birds));


%% Plot the number of "good" cells (manual and KS labelled) v. days since implant
disp(T(:,{'filename','depth','manually_sorted','probe','bad_chan','notes'}))
disp('Counting total cells')
ngood_final= [];
ngood_orig = [];

for ii = 1:height(T)
    % Find the binary directory
    root_dir = fullfile(data_dir, T.filename{ii});
    raw_dir_temp = dir(fullfile(root_dir, 'raw*'));
    raw_dir = fullfile(raw_dir_temp.folder, raw_dir_temp.name);
    ks_dir = ks_dir_fun(root_dir);

    % Latest labels saved by Phy
    [unit_ID,cluster_labels] = get_phy_cluster_labels(ks_dir);
    ngood_final(ii) = nnz(strcmp(cluster_labels,'good'));

    % Original good labels (stored in 'rez.mat')
    good_orig = getfield(load(fullfile(ks_dir,'rez.mat'),'good'),'good');
    ngood_orig(ii) = nnz(good_orig);
    
end

% Plot the number of "good" cells per session
days_since_surgery = days(T.date-T.surgery);
figure;
hold on
% cs = hsv(length(birds)+1);
hs = gobjects(length(birds),2);
for jj = 1:length(birds)
    mask =strcmp(T.bird, birds{jj});
    mask_manual = strcmp(T.manually_sorted,'yes');
    hs(jj,1) = plot(days_since_surgery(mask), ngood_orig(mask),'-o','Color',cs(jj,:),'MarkerFaceColor','w','MarkerEdgeColor',cs(jj,:));
    if any(mask_manual&mask)
    plot(days_since_surgery(mask&mask_manual), ngood_final(mask&mask_manual),'o','Color',cs(jj,:),'MarkerFaceColor',cs(jj,:),'MarkerEdgeColor',cs(jj,:));
    end
    hold on;
end
xlabel('Days since surgery'); ylabel('# "good" cells')
ylim([0 max(ngood_orig)+2])
grid on
legend(hs(:,1),birds)


%% Plot fraction E v. depth of etrode tip - good+ mua
nE = []; nI = [];
 for ii = 1:height(T)
     ks_dir = ks_dir_fun(fullfile(data_dir, T.filename{ii}));
     wvStruct = getfield(load(fullfile(ks_dir, 'waveformStruct.mat')), 'wvStruct');     
     nE(ii) = nnz(ismember(wvStruct.goodLabels,{'good','mua'}) & strcmp(wvStruct.typeLabels,'E'));
     nI(ii) = nnz(ismember(wvStruct.goodLabels,{'good','mua'}) & strcmp(wvStruct.typeLabels,'I'));
 end
 T.nE = nE(:); T.nI= nI(:);
T.ratioI = T.nI./(T.nE+T.nI);
disp(T(:,{'filename','depth','nE','nI','ratioI'}))

figure;
hs = gobjects(length(birds),1);
for ii = 1:length(birds)
    mask = strcmp(T.bird,birds{ii});
    hs(ii) = plot(T.depth(mask), T.ratioI(mask),'o','Color',...
        cs(strcmp(birds, birds{ii}),:),'LineWidth',2); hold on
%         hs(ii) = plot(T.depth(mask), T.nE(mask)+T.nI(mask),'o','Color',c(strcmp(birds, birds{ii}),:),'LineWidth',2); hold on
end
xlabel('Depth (tip, mm)')
ylabel('Fraction of cells I')
% ylabel('Total cells E+I')
grid
legend(birds)
ylim([0 .5])
