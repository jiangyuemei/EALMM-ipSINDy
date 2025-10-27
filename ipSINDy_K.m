function Xi = ipSINDy_K(Theta, dXdt, lambda_, K, n)
    % ipSINDy_K Inner-Product Sparse Identification of Nonlinear Dynamics with K-sparsity
    % 
    % Enhanced sparse regression algorithm combining thresholding and inner-product
    % driven feature selection for robust system identification under sparsity constraints
    %
    % Inputs:
    %   Theta   - Feature library matrix 
    %   dXdt    - Time derivative data or vector
    %   lambda_ - Sparsification threshold parameter  
    %   K       - Maximum sparsity level (non-zero coefficients per equation)
    %   n       - Number of state variables (system dimension)
    %
    % Output:
    %   Xi      - Sparse coefficient matrix 

    % Ensure consistent matrix formatting for dXdt
    if isvector(dXdt)
        dXdt = dXdt(:);  % Convert to column vector format
    end

    % Initialize with ordinary least squares solution
    Xi = Theta \ dXdt;  % coefficient matrix

    % Main iterative refinement loop with sparsity enforcement
    for k = 1:K
        % Apply hard thresholding based on lambda_ parameter
        smallinds = abs(Xi) < lambda_;
        Xi(smallinds) = 0;

        % Process each state variable equation independently
        for ind = 1:n
            % Identify significant feature indices for current state
            biginds = ~smallinds(:, ind);

            % Extract relevant feature subset
            Theta_big = Theta(:, biginds);
            if isempty(Theta_big)
                continue;  % Skip if no features selected
            end
            
            % Compute least squares solution on significant features
            Xi_big = Theta_big \ dXdt(:, ind);

            % Calculate residual inner products for feature ranking
            res = dXdt(:, ind);
            product = Theta_big' * res;

            % Select top-K features based on inner product magnitude
            [~, sortIdx] = sort(abs(product), 'descend');
            topK_idx = sortIdx(1:min(K, length(sortIdx)));

            % Create sparse solution preserving only top-K coefficients
            temp = zeros(size(Xi_big));
            temp(topK_idx) = Xi_big(topK_idx);

            % Update coefficient matrix with refined sparse solution
            Xi(:, ind) = 0;
            Xi(biginds, ind) = temp;
        end
    end
end
