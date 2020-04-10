clear all;

% Illustration of robust empirical Bayes confidence intervals (EBCIs)

% Reference: Armstrong, Timothy B., Michal Kolesár, and Mikkel
% Plagborg-Møller (2020), "Robust Empirical Bayes Confidence Intervals"

% See also R package "ebci": https://github.com/kolesarm/ebci

% This version: 2020-04-09
% Tested in Matlab R2019b on Windows 10 PC, 64-bit


%% Robust EBCI critical values

addpath('../ebci'); % Add folder with EBCI functions to path

% Illustration of robust EBCI critical values
disp('Critical values');
disp(cva(0, Inf, 0.05)); % m_2=0, kappa=Inf: usual normal critical value
disp(cva(4, Inf, 0.05)); % m_2=4, kappa=Inf: larger critical value that takes bias into account
disp(cva(4, 3, 0.05));   % m_2=4, kappa=3  : imposing a bound on kurtosis tightens the critical value


%% Data for small empirical example

% Load Chetty & Hendren (2018) data
dat = readtable('cz.csv', 'TreatAsEmpty', 'NA');

% For simplicity, only use 20 largest commuting zones in what follows
dat = sortrows(dat, 'pop', 'descend');
dat = dat(1:20,:);
czname = dat.czname;

% Preliminary estimates of causal effect Y_i: fixed effects estimates of
% neighborhood effect, for children with parents at the 25th percentile of
% the income distribution
Y = dat.theta25;
sigma = dat.se25; % Corresponding standard errors

% Shrink toward average outcome for permanent residents (stayers) at the
% 25th percentile of the income distribution
n = length(Y); % Sample size
X = [ones(n,1) dat.stayer25]; % Regression matrix for shrinkage
% Alternative choices for X:
% X = []; % Shrink toward zero
% X = ones(n,1); % Shrink toward grand mean

% Significance level
alpha = 0.1;

% Weights for estimating moments of effects distribution
weights = []; % empty: use equal weights


%% Parametric EBCIs

% Compute EBCIs
[thetahat, ci, w_estim, normlng, mu2, kappa, delta] = ebci(Y, X, sigma, weights, alpha, false);
% The last argument means:
% Use baseline shrinkage method that imposes moment independence (see paper)
% rather than t-statistic shrinkage

disp('Moment estimates');
disp([sqrt(mu2) kappa]);
disp('Shrinkage regression coefficients');
disp(delta');

% Compute worst-case non-coverage probability of individual parametric EBCIs
poolobj = parpool; % Open parallel computing pool
maxnoncov = nan(n,1);
parfor i=1:n % Loop over observations
    maxnoncov(i) = parametric_ebci_maxnoncov(mu2/sigma(i)^2, kappa, alpha);
end
% In this example the maximal non-coverage rates are not much larger than
% the nominal significance level of 10%, since the estimate of kappa is
% close to 3 (as in the normal distribution case) in this example

% Table with results
% thetahat: shrinkage point estimates
% ci_lo: lower endpoint of EBCIs
% ci_up: upper endpoint of EBCIs
% w_estim: shrinkage factors
% normlng: half-length of EBCIs, divided by standard errors
% maxnoncov: worst-case non-coverage probabilities
ci_lo = ci(:,1);
ci_up = ci(:,2);
ebci_param = table(czname, thetahat, ci_lo, ci_up, w_estim, normlng, maxnoncov);
disp(ebci_param);


%% Robust EBCIs, baseline shrinkage

% Call main EBCI function (uses parallel computing)
[thetahat, ci, w_estim, normlng] = ebci(Y, X, sigma, weights, alpha, false, true, true);
% Since we supply the last two arguments, we obtain robust EBCIs instead of parametric EBCIs
% The last two arguments mean (from left to right):
% - Center the EBCIs at mean-squared-error-optimal EB point estimator
% - Use estimate of kurtosis of effects distribution (if false, we only use second moment to compute critical values)

ci_lo = ci(:,1);
ci_up = ci(:,2);
ebci_robust = table(czname, thetahat, ci_lo, ci_up, w_estim, normlng);
disp(ebci_robust);

% Note that the robust EBCIs are centered at the same shrinkage estimates as the parametric EBCIs

% The robust EBCIs are only marginally wider than the parametric EBCIs,
% because the estimate of kappa is close to 3


%% Optimal robust EBCIs, baseline shrinkage

% Compute EBCIs (uses parallel computing)
[thetahat, ci, w_estim, normlng] = ebci(Y, X, sigma, weights, alpha, false, false, true);
% Penultimate argument means:
% Center the EBCIs at length-optimal EB point estimator

ci_lo = ci(:,1);
ci_up = ci(:,2);
ebci_robust_opt = table(czname, thetahat, ci_lo, ci_up, w_estim, normlng);
disp(ebci_robust_opt);

% The length-optimal robust EBCIs are only slightly shorter than the ones
% computed above, again because kappa is close to 3


%% Robust EBCIs, t-statistic shrinkage

% Compute EBCIs (does not require parallel computing due to t-stat shrinkage)
[thetahat, ci, w_estim, normlng] = ebci(Y, X, sigma, weights, alpha, true, true, true);
% Third-to-last argument means:
% Use t-statistic shrinkage, which does not impose moment independence

% Note that there's a single value of w_estim and normlng that applies to
% all observations when using t-stat shrinkage
w_estim = repmat(w_estim, n, 1);
normlng = repmat(normlng, n, 1);

ci_lo = ci(:,1);
ci_up = ci(:,2);
ebci_robust_tstat = table(czname, thetahat, ci_lo, ci_up, w_estim, normlng);
disp(ebci_robust_tstat);


% Shut down parallel computing pool
delete(poolobj);
