function [w_eb, lngth, z] = parametric_ebci(mu2_over_sigmasq, alpha)

    % Parametric Empirical Bayes confidence interval
    
    % Inputs:
    % mu2       1 x 1       signal-to-noise ratio mu_2/sigma^2
    % alpha     1 x 1       significance level
    
    % Outputs:
    % w_eb      1 x 1       shrinkage factor
    % lngth     1 x 1       half-length of CI, divided by sigma
    % z         1 x 1       critical value
    
    
    w_eb = mu2_over_sigmasq./(1+mu2_over_sigmasq); % Shrinkage factor
    z = norminv(1-alpha/2);
    lngth = sqrt(w_eb)*z; % Length

end