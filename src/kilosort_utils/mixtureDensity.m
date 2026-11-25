function D = mixtureDensity(F, gm, nu)
% Density of each row of F under each component of gm: Gaussian, or multivariate t with nu degrees of
% freedom (location gm.mu, scale gm.Sigma) when nu is given
K = gm.NumComponents; d = size(F, 2); D = zeros(size(F, 1), K);
for k = 1:K
    S = gm.Sigma(:, :, min(k, end));
    if isempty(nu)
        D(:, k) = mvnpdf(F, gm.mu(k, :), S);
    else
        Xc = F - gm.mu(k, :); delta = sum((Xc / S) .* Xc, 2);
        D(:, k) = exp(gammaln((nu + d)/2) - gammaln(nu/2) - d/2*log(nu*pi) - 0.5*log(det(S)) - (nu + d)/2*log1p(delta/nu));
    end
end
end
