function [A, B] = BDF(N, s)
    % Backward Differentiation Formula (BDF), step s ∈ [1,5]
    coeffs = {
        [1, -1],                             [1, 0];                          % s = 1
        [1, -4/3, 1/3],                      [2/3, 0, 0];                     % s = 2
        [1, -18/11, 9/11, -2/11],            [6/11, 0, 0, 0];                 % s = 3
        [1, -48/25, 36/25, -16/25, 3/25],    [12/25, 0, 0, 0, 0];             % s = 4
        [1, -300/137, 300/137, -200/137, ...
         75/137, -12/137],                  [60/137, 0, 0, 0, 0, 0];          % s = 5
    };
    a = coeffs{s, 1}; b = coeffs{s, 2};
    A = zeros(N - s + 1, N + 1);
    B = zeros(N - s + 1, N + 1);
    for i = 1:N - s + 1
        A(i, i:i+length(a)-1) = fliplr(a);
        B(i, i:i+length(b)-1) = fliplr(b);
    end
end

