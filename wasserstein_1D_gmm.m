function W = wasserstein_1D_gmm(w1, mu1, s1, w2, mu2, s2, p, N)
% WASSERSTEIN_1D_GMM  True 1D Wasserstein distance between two GMMs
%
%   W = wasserstein_1D_gmm(w1, mu1, s1, w2, mu2, s2, p, N)
%
% Inputs:
%   w1, mu1, s1 : weights, means, stds of GMM P
%   w2, mu2, s2 : weights, means, stds of GMM Q
%   p           : Wasserstein order (p = 1 or p = 2)
%   N           : number of quantile points (e.g. 500–2000)
%
% Output:
%   W           : Wasserstein-p distance
%
% Notes:
%   - Exact in 1D (up to numerical integration)
%   - No Gaussian approximation
%   - Works for arbitrary weights and variances

    arguments
        w1 (:,1) double
        mu1 (:,1) double
        s1  (:,1) double
        w2 (:,1) double
        mu2 (:,1) double
        s2  (:,1) double
        p   (1,1) double {mustBeMember(p,[1 2])}
        N   (1,1) double {mustBePositive} = 1000
    end

    % Normalize weights (safety)
    w1 = w1 / sum(w1);
    w2 = w2 / sum(w2);

    % Quantile grid (avoid 0 and 1)
    t = linspace(0,1,N+2);
    t = t(2:end-1);

    % Preallocate
    x = zeros(size(t));
    y = zeros(size(t));

    % Global search bounds
    xmin = min([mu1 - 10*s1; mu2 - 10*s2]);
    xmax = max([mu1 + 10*s1; mu2 + 10*s2]);

    % Quantile computation
    for k = 1:length(t)
        x(k) = invert_cdf(t(k), w1, mu1, s1, xmin, xmax);
        y(k) = invert_cdf(t(k), w2, mu2, s2, xmin, xmax);
    end

    % Wasserstein integral
    switch p
        case 1
            W = trapz(t, abs(x - y));
        case 2
            W = sqrt(trapz(t, (x - y).^2));
    end
end

% ---------- Helper functions ----------

function F = gmm_cdf(x, w, mu, s)
    F = 0;
    for i = 1:length(w)
        F = F + w(i) * normcdf(x, mu(i), s(i));
    end
end

function xq = invert_cdf(t, w, mu, s, xmin, xmax)
    f = @(x) gmm_cdf(x, w, mu, s) - t;
    xq = fzero(f, [xmin, xmax]);
end
