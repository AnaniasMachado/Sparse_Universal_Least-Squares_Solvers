% base_dir = fileparts(mfilename('fullpath'));

% A_filepath = fullfile(base_dir, "Experiment_Matrices", "Experiment_6", "experiment_6_matrix_m500_n250_r25_d10_idx1.mat");
expn = 7;
r_values = [25, 35];
d_values = [0.1, 0.25];
m = 100;
n = 50;
n_mtx = 5;
% output_dir = base_dir;
output_dir = "";

generate_experiment_matrices(expn, m, n, r_values, d_values, n_mtx, output_dir);