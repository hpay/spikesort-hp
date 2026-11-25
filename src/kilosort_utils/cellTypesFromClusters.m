function [types, type_p, ei_fit_p] = cellTypesFromClusters(gm, idx, P, basis)
% CELLTYPESFROMCLUSTERS  Cell type of each unit (E / I / nonsomatic / noise) from its cluster
% (gm.cluster_type; idx = 0, no features, = noise). If the model has gm.ei, units in E or I clusters
% (the somatic trough clusters; good and mua) are instead split into E and I by gm.ei: a 2-component
% mixture of t-distributions (nu = gm.ei.nu) on the per-unit values in basis named by gm.ei.features
% (trough-to-peak width in ms, asymmetry, log10 rate; from cellTypeFeatures), fit on good somatic units by
% fitCellTypeModel; component gm.ei.e_comp (the wider one) is E.
% type_p: posterior of the unit's own cluster (P from applyCellTypeModel), or for E / I units with gm.ei, of its
% E / I component. Step 7 of runSingleKilosortOpenEphys keeps units with low type_p out of good.
% ei_fit_p: for E / I units with gm.ei, how well the unit fits its own E / I component: p of its squared
% distance from the component (/ d ~ F(d, nu)), as for the cluster outliers; small = far out. NaN otherwise.

n = numel(idx); types = repmat({'noise'}, n, 1); type_p = NaN(n, 1); ei_fit_p = NaN(n, 1);
ok = idx > 0;
types(ok) = gm.cluster_type(idx(ok));
type_p(ok) = P(sub2ind(size(P), find(ok), idx(ok)));
if isfield(gm, 'ei')
    ei = ismember(types, {'E', 'I'});
    X = cell2mat(cellfun(@(f) basis.(f)(ei), gm.ei.features, 'UniformOutput', false));
    L = gm.ei.gm.ComponentProportion .* mixtureDensity(X, gm.ei.gm, gm.ei.nu);
    P_ei = L ./ sum(L, 2);
    is_e = P_ei(:, gm.ei.e_comp) >= 0.5;
    t = repmat({'I'}, nnz(ei), 1); t(is_e) = {'E'};
    types(ei) = t;
    [type_p(ei), comp] = max(P_ei, [], 2);
    d = size(X, 2); p_fit = NaN(size(comp));
    for c = 1:gm.ei.gm.NumComponents
        in_c = comp == c; Xc = X(in_c, :) - gm.ei.gm.mu(c, :);
        p_fit(in_c) = 1 - fcdf(sum((Xc / gm.ei.gm.Sigma(:, :, c)) .* Xc, 2) / d, d, gm.ei.nu);
    end
    ei_fit_p(ei) = p_fit;
end
end
