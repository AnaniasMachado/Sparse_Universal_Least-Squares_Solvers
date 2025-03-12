function A = generate_experiment_matrix(m, n, r, density, M)
    % Checks if M is given, it not define as 2
    if nargin < 5
        M = 2;
    end

    % seed = randi([1, 1000000]);

    % rng(seed);

    % Initializes vector rc
    rc = zeros(1, r);

    % Calculates values of rc
    for i = 1:r
        rc(i) = M * ((1/M) ^ (2/(i+1)));
    end

    % Generates random matrix A
    A = sprand(m, n, density, rc);
end