% A_filepath = "./Experiment_Matrices/Experiment_6/experiment_6_matrix_m500_n250_r25_d10_idx1.mat";
base_dir = fileparts(mfilename('fullpath'));

A_filepath = fullfile(base_dir, "Experiment_Matrices", "Experiment_6", "experiment_6_matrix_m500_n250_r25_d10_idx1.mat");
r = 125;
m = 500;
n = 250;
R = 1:100;
C = 1:200;
func_name = "LSFI_Det";
save_path = "Local_Search_Matrices/Experiment_6";

call_local_search_procedure(A_filepath, r, m, n, R, C, func_name, save_path);