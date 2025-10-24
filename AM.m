function [A, B] = AM(N, s)
% AM - Construct Adams-Moulton linear multistep method matrices
%
% Generates the coefficient matrices for the implicit Adams-Moulton method
% for numerical integration of ordinary differential equations.
%
% Inputs:
%   N - Number of time steps
%   s - Method order (1 to 5)
%
% Outputs:
%   A - Left-hand side coefficient matrix (N-s+1 × N+1)
%   B - Right-hand side coefficient matrix (N-s+1 × N+1)
%
% Method Orders:
%   s=1: Backward Euler (1st order)
%   s=2: Trapezoidal rule (2nd order) 
%   s=3: 3rd order Adams-Moulton
%   s=4: 4th order Adams-Moulton
%   s=5: 5th order Adams-Moulton
%
% Example:
%   [A, B] = AM(100, 3); % 3rd order Adams-Moulton for 100 steps
%
% Reference:
%   Hairer, E., & Wanner, G. (1996). Solving Ordinary Differential Equations II

    % Validate input parameters
    if s < 1 || s > 5
        error('AM:InvalidOrder', 'Method order s must be between 1 and 5');
    end
    
    if N < s
        error('AM:InsufficientSteps', ...
              'Number of steps N=%d must be >= method order s=%d', N, s);
    end

    % Adams-Moulton coefficients for orders 1 through 5
    coeffs = {
        [1, -1],                             [1/2, 1/2];                      % s = 1
        [1, -1, 0],                          [5/12, 8/12, -1/12];             % s = 2
        [1, -1, 0, 0],                       [9/24, 19/24, -5/24, 1/24];      % s = 3
        [1, -1, 0, 0, 0],                    [251/720, 646/720, -264/720, ...
                                               106/720, -19/720];             % s = 4
        [1, -1, 0, 0, 0, 0],                 [95/288, 1427/1440, -133/240, ...
                                               241/720, -173/1440, 3/160];    % s = 5
    };
    
    % Extract coefficients for specified order
    a = coeffs{s, 1};
    b = coeffs{s, 2};
    
    % Initialize coefficient matrices
    numEquations = N - s + 1;
    A = zeros(numEquations, N + 1);
    B = zeros(numEquations, N + 1);
    
    % Construct Toeplitz-like matrices for the linear multistep method
    for i = 1:numEquations
        A(i, i:i + length(a) - 1) = fliplr(a);
        B(i, i:i + length(b) - 1) = fliplr(b);
    end
    
end
