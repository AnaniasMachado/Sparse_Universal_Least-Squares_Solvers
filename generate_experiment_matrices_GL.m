function generate_experiment_matrices_GL(expn, m, r_values, d_values, n_mtx, output_dir)
    % Creates output directory if it does not exist
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % seed = randi([1, 1000000]);
    seed = 2;

    rng(seed);

    % Iterates over all values of r and d
    for i = 1:length(r_values)
        for j = 1:length(d_values)
            r = r_values(i);
            d = d_values(j);
            n = r + floor(r / 3);
            
            % Generates n_mtx matrices for each combination of r and r
            for k = 1:n_mtx
                % Generates matrix using function generate_experiment_matrix
                matrix = generate_experiment_matrix_GL(m, r, d);
                
                % Creates matrix file name
                d_decimal = round(d * 100);
                filename = sprintf('experiment_%d_matrix_m%d_n%d_r%d_d%d_idx%d.mat', expn, m, n, r, d_decimal, k);
                filepath = fullfile(output_dir, filename);
                
                % Saves matrix as a file .mat
                save(filepath, 'matrix');
                fprintf('Salvo: %s', filepath);
            end
        end
    end
end