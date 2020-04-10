function [mu2, kappa] = moment_conv(Y, sigma, weights)

    % Estimate moments of a random variable theta_i
    % Y_i = theta_i + U_i, where U ~ N(0,sigma_i^2)
    
    % Inputs:
    % Y         n x 1       data vector
    % sigma     n x 1       standard deviations
    % weights   n x 1       weight vector (if [], equal weights)
    
    % Outputs:
    % mu2       1 x 1       estimate of E[theta_i^2]
    % kappa     1 x 1       estimate of E[theta_i^4]/(E[theta_i^2]^2)

    % Determine weights
    if isempty(weights)
        weights = ones(size(Y));
    end
    weights = weights/sum(weights);
    
    % Moment estimates
    mu2 = weights'*(Y.^2-sigma.^2);
    mu4 = weights'*(Y.^4-6*(sigma.*Y).^2+3*sigma.^4);
    kappa = mu4/mu2^2;

end