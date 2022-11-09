function [idx,nlogL,P,logpdf, goodIDs] = classifyUnitsGMM(all_dir, option_only_good)

[stats, goodIDs] = getGMM_stats(all_dir, option_only_good); 
gm = getfield(load('latestGMM_cellType.mat','gm'),'gm');
[idx,nlogL,P,logpdf] = gm.cluster(table2array(stats));
