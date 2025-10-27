function [W, sigma_est, its] = AMAF(x, fx_obs, max_points, init_m_fac, max_filter_fac, expand_fac, maxits, deriv_tol, verbose)
% Computes adaptive moving average filter parameters through iterative optimization
% of bandwidth selection based on noise characteristics and signal curvature

    % Initialize algorithm parameters with default values
    if isempty(max_points)
        max_points = 10^5;
    end
    if isempty(init_m_fac)
        init_m_fac = 200;
    end
    if isempty(max_filter_fac)
        max_filter_fac = 8;
    end
    if isempty(expand_fac)
        expand_fac = 2;
    end
    if isempty(maxits)
        maxits = 100;
    end
    if isempty(deriv_tol)
        deriv_tol = 10^-6;
    end
    if isempty(verbose)
        verbose = 0;
    end

    % Estimate noise standard deviation from observed data
    sigma_est = estimate_sigma(fx_obs);
    
    % Determine subsampling factor for computational efficiency
    subsamp = max(floor(length(x) / max_points), 1);  
    
    % Create subsampled dataset
    fx_subsamp = fx_obs(1:subsamp:end, :);
    M = size(fx_subsamp, 1);
    dx_subsamp = mean(diff(x(1:subsamp:end)));
    
    % Initialize filter bandwidth parameters
    m = ceil(M / init_m_fac);
    max_filter_width = floor(M / max_filter_fac);
    
    % Initialize iteration control variables
    its = 1; 
    check = 1; 
    m = min(m, max_filter_width);

    % Main optimization loop
    while and(check > 0, its < maxits)
        % Start timer for verbose output if enabled
        if verbose
            tic
        end
        
        % Construct polynomial differentiation kernel
        [~, A] = build_poly_kernel(2, @(x) x * 0 + 1, ...
                                  min(max(floor(m * expand_fac), 3), floor((M - 1) / 2)), ...
                                  dx_subsamp, 1);
        
        % Estimate signal curvature using second derivative
        if size(fx_subsamp, 2) == 1
            d = 2 * mean(abs(conv(fx_subsamp, A(3, :), 'valid')));
        else
            d = 2 * mean(reshape(abs(conv2(A(3, :), 1, fx_subsamp, 'valid')), [], 1));
        end
        
        % Compute optimization criterion
        C = sigma_est^2 / ((d + deriv_tol)^2 * dx_subsamp^4 / 144);
        
        % Update filter bandwidth using root-finding
        mnew = min(floor((fzero(@(n) n.^5 - n.^3 - C, 1) - 1) / 2), max_filter_width);
        
        % Check convergence condition
        check = abs(m - mnew);
        m = mnew;
        its = its + 1;
        
        % Display iteration progress if verbose mode enabled
        if verbose
            disp([toc m d])
        end
    end
    
    % Scale bandwidth back to original sampling rate and construct filter
    m = m * subsamp;
    W = 1./(2 * m + 1) * ones(2 * m + 1, 1);
    
end

function sig = estimate_sigma(f)
% ESTIMATE_SIGMA - Noise standard deviation estimation via multi-scale difference analysis
%
% Estimates noise standard deviation using Richardson extrapolation principles
% applied to finite difference approximations at multiple scales

    % First-level difference computation using central differences
    Ih = (f(4:end-1) - f(2:end-3)) / 2;
    I2h = (f(5:end) - f(1:end-4)) / 4;    
    sig = sqrt(8/5) * rms(Ih - I2h);

    % Refined estimation for low-noise scenarios using higher-order differences
    if sig < 0.01
        I4h = (f(9:end) - f(1:end-8)) / 8;
        Ih_1 = 4/3 * (Ih - 1/4 * I2h);
        I2h_1 = 4/3 * (I2h(3:end-2) - 1/4 * I4h);
        sig = sqrt(576/714) * rms(Ih_1(3:end-2) - I2h_1);
    end 
    
end

function [f, A] = build_poly_kernel(deg, k, n, dx, max_dx)
% BUILD_POLY_KERNEL - Construct polynomial differentiation kernel
%
% Builds weighted polynomial basis for numerical differentiation

    % Create coordinate system and polynomial basis
    x = (-n:n)' * dx;
    X = x.^(0:deg);
    
    % Compute and normalize kernel weighting function
    K = k(x / (n * dx));
    K = K / norm(K, 1);
    
    % Compute weighted pseudoinverse for differentiation
    A = pinv(sqrt(K) .* X) .* sqrt(K)';
    
    % Construct differentiation operator matrix
    M = [diag(factorial(0:max_dx)) zeros(max_dx + 1, deg - max_dx)];
    f = M * A;
    
end