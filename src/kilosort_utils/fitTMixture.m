function gm = fitTMixture(X, K, nu)
% EM for a mixture of K multivariate t-distributions with nu degrees of freedom (full covariances),
% starting from a full-covariance Gaussian fit; clusters numbered by their mean on the first feature
[n, d] = size(X);
rng(2); g0 = fitgmdist(X, K, 'CovarianceType', 'full', 'Replicates', 20, 'RegularizationValue', 1e-3, 'Options', statset('MaxIter', 1000));
[~, order] = sort(g0.mu(:, 1)); mu = g0.mu(order, :); Sigma = g0.Sigma(:, :, order); p = g0.ComponentProportion(order);
ll_old = -Inf;
for iter = 1:2000
    L = zeros(n, K); delta = zeros(n, K);
    for k = 1:K
        Xc = X - mu(k, :); delta(:, k) = sum((Xc / Sigma(:, :, k)) .* Xc, 2);
        L(:, k) = p(k) * exp(gammaln((nu + d)/2) - gammaln(nu/2) - d/2*log(nu*pi) - 0.5*log(det(Sigma(:, :, k))) ...
            - (nu + d)/2*log1p(delta(:, k)/nu));
    end
    ll = sum(log(sum(L, 2))); R = L ./ sum(L, 2);
    U = (nu + d) ./ (nu + delta);                            % downweights units far from a component
    for k = 1:K
        w = R(:, k) .* U(:, k);
        mu(k, :) = w' * X / sum(w);
        Xc = X - mu(k, :);
        Sigma(:, :, k) = (Xc .* w)' * Xc / sum(R(:, k)) + 1e-3 * eye(d);
    end
    p = mean(R, 1);
    if ll - ll_old < 1e-7 * abs(ll); break; end
    ll_old = ll;
end
fprintf('t-mixture (%d clusters): EM converged in %d iterations\n', K, iter)
gm = gmdistribution(mu, Sigma, p);
end
