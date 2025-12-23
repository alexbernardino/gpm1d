function d = dist_kld( gm_true, gm_est )

w1 = gm_true.ComponentProportion(1);
w2 = gm_true.ComponentProportion(2);
mu1 = gm_true.mu(1);
mu2 = gm_true.mu(2);
s1 = sqrt(gm_true.Sigma(1));
s2 = sqrt(gm_true.Sigma(2));


v1 = gm_est.ComponentProportion(1);
v2 = gm_est.ComponentProportion(2);
nu1 = gm_est.mu(1);
nu2 = gm_est.mu(2);
t1 = sqrt(gm_est.Sigma(1));
t2 = sqrt(gm_est.Sigma(2));

N = 1e6;

%% --- Sample from true GMM p(x)
z = rand(N,1) < w1;
x = zeros(N,1);
x(z)  = mu1 + s1*randn(sum(z),1);
x(~z) = mu2 + s2*randn(sum(~z),1);

%% --- Log-density under p
logp = logsumexp( ...
    [log(w1) + lognormpdf(x, mu1, s1), ...
     log(w2) + lognormpdf(x, mu2, s2)], 2);

%% --- Log-density under q
logq = logsumexp( ...
    [log(v1) + lognormpdf(x, nu1, t1), ...
     log(v2) + lognormpdf(x, nu2, t2)], 2);

%% --- KL estimate
KL = mean(logp - logq);

%% --- Enforce non-negativity (numerical safeguard only)
d = max(KL, 0);

%Direct integration may not be the best approach numerically
%p = @(x) w1*normpdf(x,mu1,s1) + w2*normpdf(x,mu2,s2);
%q = @(x) v1*normpdf(x,nu1,t1) + v2*normpdf(x,nu2,t2);
%integrand = @(x) p(x) .* (log(p(x)+realmin) - log(q(x)+realmin));
%
%d = integral(integrand, -20, 20);
end

function y = lognormpdf(x, mu, sigma)
%LOGNORMPDF Log of the normal PDF, computed analytically
%
%   y = log N(x | mu, sigma^2)
%
% Numerically stable for small sigma and large |x-mu|

y = -0.5*log(2*pi) ...
    - log(sigma) ...
    - 0.5*((x - mu)./sigma).^2;
end


function y = logsumexp(a, dim)
%LOGSUMEXP Numerically stable log(sum(exp(a),dim))
%   a: array of log-values
%   dim: dimension to sum over (default 1)

if nargin < 2
    dim = 1;
end

amax = max(a,[],dim);            % max along dim
y = amax + log(sum(exp(a - amax), dim));
end
