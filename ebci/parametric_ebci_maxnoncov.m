function maxnoncov = parametric_ebci_maxnoncov(mu2_over_sigmasq, kappa, alpha, varargin)

    % Maximal non-coverage of parametric EBCI given
    % second moment mu_2 and kurtosis kappa
    
    % Inputs:
    % mu2_over_sigmasq  1 x 1       signal-to-noise ratio mu_2/sigma^2
    % kappa             1 x 1       kurtosis (set to [] or Inf if kappa=Inf)
    % alpha             1 x 1       significance level
    % varargin          struct      optional: {1} optimization options (set to [] if default options)
    
    % Outputs:
    % maxnoncov         1 x 1       maximal non-coverage probability
    
    
    % Optimization options struct
    if ~isempty(varargin)
        opt_struct = varargin{1};
    else
        opt_struct = opt_struct_default();
    end
    
    [w_eb, ~, z] = parametric_ebci(mu2_over_sigmasq, alpha);
    maxnoncov = rho((1/w_eb-1)^2*mu2_over_sigmasq, kappa, z/sqrt(w_eb), opt_struct);

end