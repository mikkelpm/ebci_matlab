function chi = cva(m2, kappa, alpha, varargin)
    
    % Robust EBCI critical value based on m_2 and possibly kappa
    
    % Inputs:
    % m2            1 x 1       second moment
    % kappa         1 x 1       kurtosis (set to [] or Inf if kappa=Inf)
    % alpha         1 x 1       significance level
    % varargin      struct      optional: {1} optimization options (set to [] if default options)
    
    % Outputs:
    % chi           1 x 1       critical value for robust EBCI

    
    if isempty(kappa)
        kappa = Inf;
    end
    
    if m2==0 || kappa==1 % Simple case
        
        chi = CVb(sqrt(m2), alpha); % Square root of non-central chi-square CV
        
    else % Otherwise, solve numerically for chi
        
        % Bounds on optimization
        lo = 0.99*CVb(sqrt(m2),alpha); % Lower bound
        up = sqrt((1+m2)/alpha); % Upper bound (from Chebyshev's inequality)
        
        % Optimization options struct
        if isempty(varargin)
            opt_struct = opt_struct_default();
        else
            opt_struct = varargin{1};
        end
        
        % First determine CV without verifying optimum
        opt_struct2 = opt_struct;
        opt_struct2.check = false;
        chi = fzero(@(cchi) rho(m2,kappa,cchi,opt_struct2)-alpha, [lo up], opt_struct.fzero);
        
        if opt_struct.check % If we should check solution, run again at putative optimum
            the_alpha = rho(m2, kappa, chi, opt_struct);
            if abs(the_alpha-alpha)>0.001
                warning('%s%f', 'The difference between non-coverage and alpha at CV is ', the_alpha-alpha);
            end
        end
    
    end

end

% Critical value from Armstrong and Kolesar
function val = CVb(B, alpha)
    idx = (B<10);
    val = B + norminv(1-alpha);
    val(idx) = sqrt(ncx2inv(1-alpha, 1, B(idx).^2));
end