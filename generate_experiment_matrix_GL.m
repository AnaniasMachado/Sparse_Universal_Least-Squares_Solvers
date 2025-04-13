function A = generate_experiment_matrix_GL(m, r, density)
    B = sprand(m, r, density);
    k = floor(r / 3);
    A = [B, B(:, 1:k) + B(:, k+1:2*k)];    
    A = full(A);
end