function A = generate_experiment_matrix(m, n, r, density)
    % Initializes vector rc
    rc = zeros(1, r);

    % Initializes constant M
    M = 2;

    % Initialites constant rho
    rho=(1/M)^(2/(r+1));

    % Specifies entries of vector rc
    for i=1:r
        rc(i)=M*rho^i;
    end

    % Generates random matrix A
    A = sprand(m, n, density, rc);
end