function [idx,labels] = GMM_apply(gm, all_dir, option_only_good, plot_on)
%%
id = 'stats:gmdistribution:cluster:MissingData'; warning('off',id)
n_clusters = gm.NumComponents;
[stats, goodIDs, mxWF, tWF] = GMM_stats(all_dir, option_only_good); 
[idx, nlogL, P, logpdf, d2] = gm.cluster(table2array(stats));
warning('on',id)
% P: posterior probabilities of each Gaussian mixture component in gm given each observation in X
% logpdf: logarithm of the estimated probability density function (pdf) evaluated at each observation in X
% d2: squared Mahalanobis distance of each observation in X to each Gaussian mixture component in gm

% Categorize cells where 1. max is greater than min or 2. stats could not be
% calculated as unknown
mask1 =  (max(mxWF,[],2)>abs(min(mxWF,[],2))) | isnan(idx);
idx(mask1) = 0;

% Categorize cells where 3. trough is not at expected time (>2 sigma) as unknown
[~, mins] = min(mxWF,[],2);
tmins = tWF(mins);
mean_t = mean(tmins(idx>0));
std_t = std(tmins(idx>0));
% mask_tmin = abs(tmins)>0.25*1e-3; % (s)
mask_tmin = abs((tmins - mean_t)/std_t)>3;
idx(mask_tmin) = 0;

% Apply labels
labels = cell(size(idx));
labels(idx==0) = {'unknown'};
labels(idx==1) = {'E'};
labels(idx==2) = {'I'};
inds_valid = find(idx>0);

% Extract relevant post. prob and Mahalanobis distances, take log to plot
P_match = NaN(size(idx));
d2_match = NaN(size(idx));

P_match(inds_valid) = P(sub2ind(size(P),inds_valid, idx(inds_valid)));
P_plot = log(P_match);
d2_match(inds_valid) = d2(sub2ind(size(d2),inds_valid, idx(inds_valid)));
d2_plot = log(d2_match);
if n_clusters==2
  % Apply posterior probability criteria
    min_P = 0.7;
    mask = P_match<min_P & idx>0;
    idx(mask) = 3;
    labels(mask) = {'intermediate'};
    
    % Apply distance criteria
    max_d2a = 10;
    mask = d2_match>max_d2a & idx>0;
    idx(mask) = 4;
    labels(mask) = {'moderate_outlier'};
    
    
    max_d2b = 20;
    mask = d2_match>max_d2b & idx>0;
    idx(mask) = 5;
    labels(mask) = {'extreme_outlier'};
  
end
% function [idx,nlogL,P,logpdf,d2] = applyGMM(gm, stats, mxWF, tWF, plot_on)
fprintf('E %i, I %i, intermediate %i, moderate_outlier %i, extreme_outlier %i, unclassified %i\n', ...
    nnz(idx==1), nnz(idx==2), nnz(idx==3), nnz(idx==4), nnz(idx==5), nnz(idx==0))
% disp(mean(stats{idx==1,:}))
% disp(mean(stats{idx==2,:}))

if plot_on
    %% Plot results
    fs = 1/mean(diff(tWF));
    figure; 
%     set(gcf,'Pos',[ 4657        -119        1920         955])
    colormap(cool)

    clf
    hs = subplot(2,4,1);
    scatter3( stats.duration/fs*1e3, stats.asymmetry,stats.rate_log10, 15,stats.pk_deriv, 'filled');   axis vis3d
    xlabel('Duration tr-pk(ms)'); ylabel('Asymmetry'); zlabel('log10(rate)')
    c=colorbar('north'); c.Label.String = 'Ratio peak:trough deriv.';
    
    hs(2) = subplot(2,4,2);
    scatter3(stats.duration/fs*1e3, stats.asymmetry,stats.rate_log10,15,idx,'filled'); axis vis3d
    xlabel('Duration tr-pk(ms)'); ylabel('Asymmetry'); zlabel('log10(rate)')
    if ischar(all_dir); title(all_dir,'Interp','none'); end
    c=colorbar('north');
    c.Ticks = 0:5;
    c.TickLabels = {'unknown','E','I','intermediate','moderate_outlier','extreme_outlier'};
    c.TickLabelInterpreter = 'none';
    
    hs(3) = subplot(2,4,3);
    scatter3(stats.duration/fs*1e3, stats.asymmetry,stats.rate_log10,15,d2_plot,'filled'); axis vis3d
    xlabel('Duration tr-pk(ms)'); ylabel('Asymmetry'); zlabel('log10(rate)')
    c=colorbar('north');
    c.Label.String = 'log squared Mahalanobis distance';
      
    hs(4) = subplot(2,4,4);
    scatter3(stats.duration/fs*1e3, stats.asymmetry,stats.rate_log10,15,P_plot,'filled'); axis vis3d
    xlabel('Duration tr-pk(ms)'); ylabel('Asymmetry'); zlabel('log10(rate)')
    c=colorbar('north');
    c.Label.String = 'log posterior P';
    %     linkaxes(hs);
    Link = linkprop(hs,{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
    setappdata(gcf, 'StoreTheLink', Link);
     axis(hs,'square')
    xlim([0  1]); ylim([-1 1]); zlim([-2.3 1.6]); set(gca,'XTick',0:.2:1)
    
    % Plot waveforms
    hs = subplot(2,4,5); 
    if any(idx==1); plot(tWF,zscore(mxWF(idx==1, :)'),'k');   hold on,plot(tWF,mean(zscore(mxWF(idx==1,:)'),2),'r','LineWidth',2); end
    axis square; title('Excitatory'); ylabel('Z-scored amplitude')    
    
    hs(2) = subplot(2,4,6); 
    if any(idx==2); plot(tWF,zscore(mxWF(idx==2, :)'),'k');  hold on,plot(tWF,mean(zscore(mxWF(idx==2,:)'),2),'b','LineWidth',2); end
    axis square; title('Inhibitory')
    
    hs(3) = subplot(2,4,7);
    if any(idx==0); plot(tWF,zscore(mxWF(idx==0, :)'),'Color','k');  hold on; end
    if any(idx==3); plot(tWF,zscore(mxWF(idx==3, :)'),'Color',[ 0.6000    0.4000    1.0000]); hold on; end
        axis square; title('Intermediate/unknown')

        hs(4) = subplot(2,4,8);
    if any(idx==4); plot(tWF,zscore(mxWF(idx==4, :)'),'Color',[0.8 .2 1]); hold on; end
    if any(idx==5); plot(tWF,zscore(mxWF(idx==5, :)'),'Color','k');  end
%     legend('intermediate', 'moderate_outlier','extreme_outlier','unknown')
    axis square; title('moderate outlier/extreme outlier')
    linkaxes(hs); xlim(tWF([1 end])); 
    ylim([-6 3])
    grid(hs,'on')
    drawnow
end

