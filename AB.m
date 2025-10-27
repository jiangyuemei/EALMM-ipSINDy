function [A, B] = AB(N, s)
    % Adams-Bashforth method (explicit), step s ∈ [1,5]
    coeffs = {
        [1, -1],                             [0, 1];                          % s = 1
        [1, -1, 0],                          [0, 3/2, -1/2];                  % s = 2
        [1, -1, 0, 0],                       [0, 23/12, -16/12, 5/12];        % s = 3
        [1, -1, 0, 0, 0],                    [0, 55/24, -59/24, 37/24, -9/24];% s = 4
        [1, -1, 0, 0, 0, 0],                 [0, 1901/720, -2774/720, ...
                                               2616/720, -1274/720, 251/720]  % s = 5
    };
    a = coeffs{s, 1}; b = coeffs{s, 2};
    A = zeros(N - s + 1, N + 1);
    B = zeros(N - s + 1, N + 1);
    for i = 1:N - s + 1
        A(i, i:i+length(a)-1) = fliplr(a);
        B(i, i:i+length(b)-1) = fliplr(b);
    end
end

