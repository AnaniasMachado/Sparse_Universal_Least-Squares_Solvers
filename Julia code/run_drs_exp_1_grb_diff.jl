using CSV
using DataFrames
using MAT
using Statistics

include("utility.jl")

exp = "1"
grb_exp = "10"

matrices_folder = "./Experiment_Matrices/DRS_Experiment_$(exp)"
mat_files_exp = "1"

methods = ["ADMM", "ADMMe", "DRS", "DRS_Boyd", "DRS_FP"]

problems = ["P123"]
problem = problems[1]

gurobi_solutions_folder = "./Solutions/Experiment_$(grb_exp)/Gurobi_Cal"

results_folder = "results/DRS_Experiment_$exp"

m_values = [60+i*20 for i in 0:27]
d_values = [100]

df = DataFrame()

for method in methods
    solutions_folder = "./Solutions/DRS_Experiment_$(exp)/$(method)"

    H_div_AMP_norm_0_diff_list = []
    H_div_AMP_norm_1_diff_list = []
    for m in m_values
        n = Int(m / 2)
        r = Int(m / 4)
        for d in d_values
            for idx in 1:5
                try
                    sol_name = "Experiment_$(exp)_$(problem)_m_$(m)_n_$(n)_idx_5"
                    sol_path = joinpath(solutions_folder, sol_name)
                    sol_data = matread(sol_path)

                    grb_name = "Experiment_$(grb_exp)_P123_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                    grb_path = joinpath(gurobi_solutions_folder, grb_name)
                    grb_data = matread(grb_path)
                catch
                    println("Caught error. m = $(m)")
                    break
                end
                mat_name = "experiment_$(mat_files_exp)_matrix_m$(m)_n$(n)_r$(r)_d$(d)_idx$(idx).mat"
                mat_path = joinpath(matrices_folder, mat_name)
                mat_data = matread(mat_path)
                A = mat_data["matrix"]
                A = Matrix(A)
                AMP = pinv(A)

                grb_name = "Experiment_$(grb_exp)_P123_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                grb_path = joinpath(gurobi_solutions_folder, grb_name)
                grb_data = matread(grb_path)
                H_GRB = grb_data["H"]
                H_GRB = Matrix(H_GRB)

                H_GRB_norm_0 = matrix_norm_0(H_GRB)
                H_GRB_norm_1 = norm(H_GRB, 1)

                H_div_AMP_norm_0_GRB = H_GRB_norm_0 / matrix_norm_0(AMP)
                H_div_AMP_norm_1_GRB = H_GRB_norm_1 / norm(AMP, 1)

                sol_name = "Experiment_$(exp)_$(problem)_m_$(m)_n_$(n)_idx_$(idx)"
                sol_path = joinpath(solutions_folder, sol_name)
                sol_data = matread(sol_path)
                H = sol_data["H"]
                H = Matrix(H)
                time = sol_data["time"]

                H_norm_0 = matrix_norm_0(H)
                H_norm_1 = norm(H, 1)

                H_div_AMP_norm_0 = H_norm_0 / matrix_norm_0(AMP)
                H_div_AMP_norm_1 = H_norm_1 / norm(AMP, 1)

                push!(H_div_AMP_norm_0_diff_list, H_div_AMP_norm_0 - H_div_AMP_norm_0_GRB)
                push!(H_div_AMP_norm_1_diff_list, H_div_AMP_norm_1 - H_div_AMP_norm_1_GRB)
            end
        end
    end

    H_div_AMP_norm_0_min = minimum(H_div_AMP_norm_0_diff_list)
    H_div_AMP_norm_0_mean = mean(H_div_AMP_norm_0_diff_list)
    H_div_AMP_norm_0_max = maximum(H_div_AMP_norm_0_diff_list)

    H_div_AMP_norm_1_min = minimum(H_div_AMP_norm_1_diff_list)
    H_div_AMP_norm_1_mean = mean(H_div_AMP_norm_1_diff_list)
    H_div_AMP_norm_1_max = maximum(H_div_AMP_norm_1_diff_list)

    result = DataFrame(
        method = [method],
        H_div_AMP_norm_0_min = [H_div_AMP_norm_0_min],
        H_div_AMP_norm_0_mean = [H_div_AMP_norm_0_mean],
        H_div_AMP_norm_0_max = [H_div_AMP_norm_0_max],
        H_div_AMP_norm_1_min = [H_div_AMP_norm_1_min],
        H_div_AMP_norm_1_mean = [H_div_AMP_norm_1_mean],
        H_div_AMP_norm_1_max = [H_div_AMP_norm_1_max]
    )

    append!(df, result)
end

results_filename = "results_$(exp)_$(problem)_grb_diff.csv"
results_filepath = joinpath(results_folder, results_filename)
CSV.write(results_filepath, df)