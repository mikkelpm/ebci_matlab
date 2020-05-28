function [mu2, kappa] = moment_conv(Y, sigma, weights, fs_correction)

    % Estimate moments of a random variable theta_i
    % Y_i = theta_i + U_i, where U ~ N(0,sigma_i^2)
    
    % Inputs:
    % Y             n x 1       data vector
    % sigma         n x 1       standard deviations
    % weights       n x 1       weight vector (if [], equal weights)
    % fs_correction char        finite-sample moment correction, either 'none', 'PMT', or 'FPLIB'
    
    % Outputs:
    % mu2       1 x 1       estimate of E[theta_i^2]
    % kappa     1 x 1       estimate of E[theta_i^4]/(E[theta_i^2]^2)

    % Determine weights
    if isempty(weights)
        weights = ones(size(Y));
    end
    weights = weights/sum(weights); % Normalize to sum to 1
    
    % Preliminary moment estimates
    W2 = Y.^2-sigma.^2;
    W4 = Y.^4-6*(sigma.*Y).^2+3*sigma.^4;
    mu2_t = weights'*W2;
    mu4_t = weights'*W4;
    
    % Finite-sample moment corrections
    wgtV = @(Z) ((weights.^2)'*(Z.^2-(weights'*Z)^2))/(1-(weights'*weights));
    tmean = @(m,V) m + sqrt(V)*normpdf(m/sqrt(V))/normcdf(m/sqrt(V)); % Mean of truncated normal
    
    switch upper(fs_correction)
        
        case 'NONE'
            mu2 = max(mu2_t, 0);
            kappa = max(mu4_t/mu2^2, 1);
        case 'PMT'
            mu2 = max(mu2_t, 2*((weights.^2)'*(sigma.^4))/(weights'*(sigma.^2)));
            kappa = max(mu4_t/mu2^2, 1 + 32*((weights.^2)'*(sigma.^8))/(mu2^2*(weights'*(sigma.^4))));
        case 'FPLIB'
            mu2 = tmean(mu2_t, wgtV(W2));
            kappa = 1 + tmean(mu4_t-mu2_t^2, wgtV(W4-2*mu2*W2))/mu2^2;
            
    end
    
end
