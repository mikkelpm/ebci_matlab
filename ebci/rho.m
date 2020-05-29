function [maxnoncov, lf_t, lf_p] = rho(m2, kappa, chi, varargin)

    % Maximal non-coverage of EBCI based on second moment m_2,
    % kurtosis kappa, and critical value chi
    
    % Inputs:
    % m2            1 x 1       second moment
    % kappa         1 x 1       kurtosis (set to [] or Inf if kappa=Inf)
    % chi           1 x 1       critical value
    % varargin      struct      optional: {1} optimization options (set to [] if default options)
    
    % Outputs:
    % maxnoncov     1 x 1       maximal non-coverage probability
    % lf_t          1 x q       support points for least-favorable (LF) distribution for t=b^2        
    % lf_p          1 x q       probability masses for LF distribution for t=b^2

    
    if isempty(kappa)
        kappa = Inf;
    end
    
    % If distribution is a point mass
    if kappa==1
        maxnoncov = r(m2, chi);
        lf_t = m2;
        lf_p = 1;
        return;
    end
    
    % Optimization options struct
    if isempty(varargin)
        opt_struct = opt_struct_default();
    else
        opt_struct = varargin{1};
    end
    
    % Otherwise, first solve problem with only second moment
    [maxnoncov, t0, ip, lf_t, lf_p] = rho0(m2, chi, opt_struct);
    
    % If the above solution already satisfies kurtosis bound, we're done
    if lf_p*(lf_t.^2)' <= kappa*m2^2
        return;
    end
    
    % Otherwise, solve the optimization that imposes kappa bound
    
    % First determine where delta(x, x_0) is maximized
    [~, tbar] = lam(0, chi, t0, ip, opt_struct);
    
    % Optimize separately below and above bar{t}, since there are
    % typically multiple local minima
    obj = @(x0) r(x0, chi) + r1(x0, chi)*(m2-x0) ...
                + lammax(x0, chi, t0, ip, tbar, opt_struct) * (kappa*m2^2-2*x0*m2+x0^2);
    
    if tbar>0
        [x0opt_below,maxnoncov_below] = fminbnd(obj, 0, tbar, opt_struct.fminbnd); % Solution below bar{t}
    else
        x0opt_below = 0;
        maxnoncov_below = obj(0);
    end
    [x0opt, maxnoncov] = fminbnd(obj, tbar, t0, opt_struct.fminbnd); % Solution above bar{t}
    
    if maxnoncov_below<maxnoncov
        x0opt = x0opt_below;
        maxnoncov = maxnoncov_below;
    end
    
    % Least favorable (LF) distribution
    [~, xopt] = lam(x0opt, chi, t0, ip, opt_struct);
    lf_t = sort([x0opt xopt]); % Support points
    p = (m2-lf_t(2))/(lf_t(1)-lf_t(2));
    lf_p = [p 1-p]; % Probability masses
    
    % Non-coverage rate, m_2, and kappa at LF solution
    primal_maxnoncov = r(lf_t, chi)*lf_p';
    primal_m2 = lf_t*lf_p';
    primal_kappa = (lf_t.^2)*lf_p'/primal_m2^2;
    
    % If LF solution is close to dual, no need to check linear program.
    % Otherwise, double-check we found the optimum by solving the primal
    if max([maxnoncov m2 kappa] - [primal_maxnoncov primal_m2 primal_kappa]) > 1e-4 ...
       && (isfield(opt_struct, 'check') && opt_struct.check)
       
        % Add m_2 here for cases where it's very small, so we can satisfy
        % the constraint.
        xs = unique([m2 lf_t linspace(0, t0, opt_struct.numgrid)], 'sorted');
        numgrid = length(xs);
        [~, mopt, exitflag] = linprog(-r(xs, chi)', ... % Objective
                            xs.^2, kappa*m2^2, ... % Inequality constraint
                            [ones(1,numgrid); xs], [1; m2], ... % Equality constraints
                            zeros(1,numgrid), [], opt_struct.linprog);
        opt = -mopt;
        
        if exitflag~=1 || abs(opt-maxnoncov)>1e-4
            warning('%s%f\n%s%f\n%s\n%s%f%s%f%s%f', ...
                    'Linear program finds non-coverage ', opt, ...
                    'Direct approach finds non-coverage ', maxnoncov, ...
                    'Difference>0.001. This happened for', ...
                    'chi=', chi, ', mu_2=', m2, ', kappa=', kappa);
        end
        
    end

end

function [maxnoncov, t0, ip, lf_t, lf_p] = rho0(t, chi, opt_struct)

    % Maximal non-coverage based only on second moment

    [t0, ip] = rt0(chi, opt_struct);
    
    if t<t0 % Least favorable distribution has two mass points
        lf_t = [0 t0];
        lf_p = [1-t/t0 t/t0];
    else % Least favorable distribution has one mass point
        lf_t = t;
        lf_p = 1;
    end
    maxnoncov = lf_p*r(lf_t,chi)';

end

function [val, xmax] = lam(x0, chi, t0, ip, opt_struct)

    % Maximize delta(x, x_0, chi) over x
    
    % Check derivatives at 0, inflection point, t0, and x0. If we're above
    % inflection point, then maximum is below it, and it's at zero if
    % derivative at zero is negative. Otherwise between 0 and inflection
    % point.
    xs = sort([t0 ip]);
    if x0 >= xs(1)
        xs = [0 xs(1)];
    else
        xs = unique([0 x0 xs]);
    end
    
    vals = delta(xs, x0, chi);
    ders = delta1(xs, x0, chi);
    
    % Default return value: optimum at 0
    val = vals(1);
    xmax = 0;
    
    % Expect delta has single maximum, so first increasing, then
    % decreasing, up to numerical tolerance
    if all(ders<=0) && vals(1)==max(vals)
        % Maximum at 0
        return;
    elseif all(diff(ders>=0)<=0) && ders(end)<=0
        % Function first increasing, then decreasing
        the_ind = max(find(ders<0,1),2); % In case all derivatives are negative
        if isempty(the_ind)
            the_ind = length(xs);
        end
        the_start = xs(the_ind-1);
        the_end = xs(the_ind);
    elseif (min(abs(ders)) < 1e-6)
        % Determine interval based on value of delta,
        % numerical accuracy of delta1 only 7e-6
        [~,the_ind_max] = max(vals);
        the_start = xs(max(the_ind_max-1, 1));
        the_end = xs(min(the_ind_max+1, length(xs)));
    else
        error('%s%f%s%f%s', 'There are multiple local optima in the function delta(x, x0=', x0, ', chi=', chi, ').');
    end
    
    % Numerical optimization
    [the_xmax, mdelta] = fminbnd(@(x) -delta(x, x0, chi), the_start, the_end, opt_struct.fminbnd);
    
    % Check optimum at 0 not substantially higher, we could miss it due to numerical accuracy issues
    if -mdelta > vals(1)
        val = -mdelta;
        xmax = the_xmax;
    elseif -mdelta < vals(1)-1e-9
        warning('%s%f%s%f%s', 'Optimum may be wrong for lam(x0=', x0, ', chi=', chi, ').');
    end

end

function val = lammax(x0, chi, t0, ip, tbar, opt_struct)
    if x0 >= tbar
        val = delta(0, x0, chi);
    else
        val = lam(x0, chi, t0, ip, opt_struct);
    end
end

function [t0, ip] = rt0(chi, opt_struct)

    % Find t_0 and inflection point, called t_1 in the paper
    
    if chi<sqrt(3)
        t0 = 0;
        ip = 0;
    else
        % Find inflection point
        % Avoid issues when chi is numerically very large
        if abs(r2(chi^2-3/2, chi))<1e-12 || (chi^2-3)==chi^2
            ip = chi^2-3/2;
        else
            ip = fzero(@(tt) r2(tt,chi), [chi^2-3, chi^2], opt_struct.fzero);
        end
        
        % Find t_0
        f0 = @(tt) r(tt,chi)-tt.*r1(tt,chi)-r(0,chi);
        % Make sure upper endpoint of interval is positive; it always is for
        % chi< 100,000, so we should never enter the while loop
        lo = ip;
        up = 2*chi^2;
        while f0(up) < 0
            lo = up;
            up = 2*up;
        end
        t0 = lo;
        if f0(lo) < 0
            t0 = fzero(f0, [lo up], opt_struct.fzero);    
        elseif f0(lo) > 1e-12
            warning('%s%f', 'Failed to solve for t0 using rt0 at chi=', chi);
        end
    end
    
end

function val = r(t, chi)
    % Function called r_0 in paper
    idx = (sqrt(t)-chi<=5);
    val = ones(size(t));
    val(idx) = normcdf(-sqrt(t(idx))-chi) + normcdf(sqrt(t(idx))-chi);
end

function val = r1(t, chi)
    % First derivative of r
    % Apply l'Hopital's rule. This gives maximum absolute error 1e-9.
    idx = (t>=1e-8);
    val = repmat(chi*normpdf(chi), size(t));
    val(idx) = (normpdf(sqrt(t(idx))-chi)-normpdf(sqrt(t(idx))+chi))./(2*sqrt(t(idx)));
end

function val = r2(t, chi)
    % Second derivative of r
    % Apply L'Hopital's rule 3x. This gives maximum abs error about 8e-8.
    idx = (t>=2e-6);
    val = repmat(normpdf(chi)*chi*(chi^2-3)/6, size(t));
    val(idx) = (normpdf(sqrt(t(idx))+chi).*(chi*sqrt(t(idx))+t(idx)+1) ...
                + normpdf(sqrt(t(idx))-chi).*(chi*sqrt(t(idx))-t(idx)-1)) ...
                ./ (4*t(idx).^(3/2));
end

function val = r3(t, chi)
    % Third derivative of r
    % Apply L'Hopital's rule 5x. This gives maximum abs error about 3e-6.
    idx = (t>=2e-4);
    val = repmat(normpdf(chi)*(chi^5-10*chi^3+15*chi)/60, size(t));
    val(idx) = (normpdf(chi-sqrt(t(idx))).*(t(idx).^2-2*chi*t(idx).^(3/2)+(2+chi^2)*t(idx)-3*chi*sqrt(t(idx))+3) ...
                - normpdf(chi+sqrt(t(idx))).*(t(idx).^2+2*chi*t(idx).^(3/2)+(2+chi^2)*t(idx)+3*chi*sqrt(t(idx))+3)) ...
                ./ (8*t(idx).^(5/2));
end

function val = delta(x, x0, chi)
    % Apply L'Hopital's rule 2x. This gives maximum abs error about 1e-7.
    idx = (abs(x-x0)>=1e-4);
    val = repmat(r2(x0, chi)/2, size(x));
    val(idx) = (r(x(idx), chi) - r(x0, chi) - r1(x0, chi).*(x(idx)-x0)) ./ (x(idx)-x0).^2;
end

function val = delta1(x, x0, chi)
    % Apply L'Hospital's rule 2x. This gives maximum abs error about 7e-6.
    idx = (abs(x-x0)>=1e-3);
    val = repmat(r3(x0, chi)/6, size(x));
    val(idx) = ((r1(x(idx), chi)+r1(x0, chi))-2*(r(x(idx), chi)-r(x0, chi))./(x(idx)-x0))./(x(idx)-x0).^2;
end