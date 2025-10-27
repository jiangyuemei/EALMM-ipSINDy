function [A, B] = AM(N, s)
    % Adams-Moulton method (implicit), step s ∈ [1,5]
    coeffs = {
        [1, -1],                             [1/2, 1/2];                      % s = 1
        [1, -1, 0],                          [5/12, 8/12, -1/12];             % s = 2
        [1, -1, 0, 0],                       [9/24, 19/24, -5/24, 1/24];      % s = 3
        [1, -1, 0, 0, 0],                    [251/720, 646/720, -264/720, ...
                                               106/720, -19/720];             % s = 4
        [1, -1, 0, 0, 0, 0],                 [95/288, 1427/1440, -133/240, ...
                                               241/720, -173/1440, 3/160];    % s = 5
    };
    a = coeffs{s, 1}; b = coeffs{s, 2};
    A = zeros(N - s + 1, N + 1);
    B = zeros(N - s + 1, N + 1);
    for i = 1:N - s + 1
        A(i, i:i+length(a)-1) = fliplr(a);
        B(i, i:i+length(b)-1) = fliplr(b);
    end
end

