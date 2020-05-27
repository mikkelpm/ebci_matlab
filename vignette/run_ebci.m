clear all;

% Illustration of robust empirical Bayes confidence intervals (EBCIs)

% Reference: Armstrong, Timothy B., Michal Kolesár, and Mikkel
% Plagborg-Møller (2020), "Robust Empirical Bayes Confidence Intervals"
% https://arxiv.org/abs/2004.03448

% See also R package "ebci": https://github.com/kolesarm/ebci

% Tested in Matlab R2019b on Windows 10 PC, 64-bit


%% Robust EBCI critical values

addpath('../ebci'); % Add folder with EBCI functions to path

% Illustration of robust EBCI critical values
disp('** Critical values **');
disp('m_2=0, kappa=Inf');
disp(cva(0, Inf, 0.05)); % Usual normal critical value
disp('m_2=4, kappa=Inf');
disp(cva(4, Inf, 0.05)); % Larger critical value that takes bias into account
disp('m_2=4, kappa=3');
disp(cva(4, 3, 0.05));   % Imposing a bound on kurtosis tightens the critical value


%% Data for small empirical example

disp('** Empirical example **');

% Load Chetty & Hendren (2018) data
dat = readtable('cz.csv', 'TreatAsEmpty', 'NA');
dat = dat(~isnan(dat.theta25),:); % Drop missing
czname = dat.czname; % Commuting zone
state = dat.state; % State

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

% Precision weights for regression and moment estimates
weights = 1./sigma.^2;

% Significance level
alpha = 0.1;

poolobj = parpool; % Open parallel computing pool


%% Robust EBCIs, baseline shrinkage

% Compute robust EBCIs (uses parallel computing)
[thetahat, ci, w_estim, normlng, mu2, kappa, delta] = ebci(Y, X, sigma, alpha, 'weights', weights);

disp('Moment estimates [mu_2 kappa]');
disp([mu2 kappa]);
disp('Shrinkage regression coefficients [intercept stayer25]');
disp(delta');

% Table with results
% thetahat: shrinkage point estimates
% ci_lo: lower endpoint of EBCIs
% ci_up: upper endpoint of EBCIs
% w_estim: shrinkage factors
% normlng: half-length of EBCIs, divided by standard errors
ci_lo = ci(:,1);
ci_up = ci(:,2);
ebci_robust = table(czname, state, thetahat, ci_lo, ci_up, w_estim, normlng);

% Display results for California
disp('Robust EBCIs: California');
disp(ebci_robust(strcmp(ebci_robust.state, 'CA'),:));

disp('Average length of unshrunk CIs relative to robust EBCIs');
disp((2*norminv(1-alpha/2)*mean(sigma))/mean(ebci_robust.ci_up-ebci_robust.ci_lo));


%% Parametric EBCIs (Morris, 1983)

% Compute parametric EBCIs
[thetahat, ci, w_estim, normlng] = ebci(Y, X, sigma, alpha, 'param', true, 'weights', weights);
% The argument 'param'=true asks for parametric EBCIs
% Note that the parametric EBCIs are centered at the same shrinkage estimates as the robust EBCIs

% Compute worst-case non-coverage probability of individual parametric EBCIs
maxnoncov = nan(n,1);
parfor i=1:n % Loop over observations
    maxnoncov(i) = parametric_ebci_maxnoncov(mu2/sigma(i)^2, kappa, alpha);
end

ci_lo = ci(:,1);
ci_up = ci(:,2);
ebci_param = table(czname, state, thetahat, ci_lo, ci_up, w_estim, normlng, maxnoncov);

% Display results for California
disp('Parametric EBCIs: California');
disp(ebci_param(strcmp(ebci_param.state, 'CA'),:));

disp('Average length of parametric EBCIs relative to robust EBCIs');
disp(mean(ebci_param.ci_up-ebci_param.ci_lo)/mean(ebci_robust.ci_up-ebci_robust.ci_lo));

disp('Average worst-case non-coverage probability of parametric EBCIs');
disp(mean(maxnoncov));

% The parametric EBCIs are shorter but less robust: They could potentially
% have Empirical Bayes non-coverage probabilities well above alpha=10%


% Shut down parallel computing pool
delete(poolobj);
