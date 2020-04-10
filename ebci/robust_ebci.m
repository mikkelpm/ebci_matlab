function [w_estim, normlng] = robust_ebci(w, mu2_over_sigmasq, kappa, alpha, varargin)

    % Robust Empirical Bayes confidence interval
    
    % Optionally computes EBCI-length-optimal shrinkage factor w_opt
    
    % Inputs:
    % w                 1 x 1       shrinkage weight (set to [] if we should compute length-optimal w_opt)
    % mu2_over_sigmasq  1 x 1       signal-to-noise ratio mu_2/sigma^2
    % kappa             1 x 1       kurtosis (set to [] or Inf if kappa=Inf)
    % alpha             1 x 1       significance level
    % varargin          struct      optional: {1} optimization options (set to [] if default options)
    
    % Outputs:
    % w_estim           1 x 1       equal to w if not empty, otherwise equal to w_opt
    % normlng           1 x 1       half-length of robust EBCI, divided by sigma
    
    
    % Optimization options struct
    if isempty(varargin)
        opt_struct = opt_struct_default();
    else
        opt_struct = varargin{1};
    end

    % Function returning normalized half-length of robust EBCI
    length_fct = @(ww) ww*cva((1/ww-1)^2*mu2_over_sigmasq, kappa, alpha, opt_struct);
    
    if isempty(w) % If w is not pre-specified, find w=w_opt that minimizes CI length
        [w_estim, normlng] = fminbnd(length_fct, mu2_over_sigmasq/(1+mu2_over_sigmasq), 1, opt_struct.fminbnd);
    else % Otherwise, use critical value corresponding to pre-specified w
        w_estim = w;
        normlng = length_fct(w);
    end

end