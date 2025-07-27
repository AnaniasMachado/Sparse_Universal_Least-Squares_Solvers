using CSV
using DataFrames
using MAT
using Statistics

include("utility.jl")

exp = "3"
grb_exp = "11"

matrices_folder = "./Experiment_Matrices/DRS_Experiment_$(exp)"
mat_files_exp = "5"

methods = ["Gurobi", "ADMM", "ADMMe", "DRS", "DRS_Boyd", "DRS_FP"]

problems = ["P134"]
problem = problems[1]

gurobi_solutions_folder = "./Solutions/Experiment_$(grb_exp)/Gurobi"

results_folder = "results/DRS_Experiment_$exp"

m_values = [60+i*20 for i in 0:27]
d_values = [100]

df = DataFrame()

for method in methods
    solutions_folder = "./Solutions/DRS_Experiment_$(exp)/$(method)"

    H_div_bound_norm_0_list = []
    H_div_AMP_norm_0_list = []
    H_div_AMP_norm_1_list = []
    H_rank_ratio_list = []
    m_max = 600
    for m in m_values
        n = m
        r = Int(m / 4)
        for d in d_values
            for idx in 1:5
                try
                    if method != "Gurobi"
                        sol_name = "Experiment_$(exp)_$(problem)_m_$(m)_n_$(n)_idx_5"
                        sol_path = joinpath(solutions_folder, sol_name)
                        sol_data = matread(sol_path)
                    else
                        grb_name = "Experiment_$(grb_exp)_P13R_P14R_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                        grb_path = joinpath(gurobi_solutions_folder, grb_name)
                        grb_data = matread(grb_path)
                    end
                catch
                    println("Caught error. m = $(m)")
                    m_max = min(m, m_max)
                    break
                end
                mat_name = "experiment_$(mat_files_exp)_matrix_m$(m)_n$(n)_r$(r)_d$(d)_idx$(idx).mat"
                mat_path = joinpath(matrices_folder, mat_name)
                mat_data = matread(mat_path)
                A = mat_data["matrix"]
                A = Matrix(A)
                AMP = pinv(A)

                H = 0
                if method == "Gurobi"
                    grb_name = "Experiment_$(grb_exp)_P13R_P14R_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                    grb_path = joinpath(gurobi_solutions_folder, grb_name)
                    grb_data = matread(grb_path)
                    H = grb_data["H"]
                    H = Matrix(H)
                else
                    sol_name = "Experiment_$(exp)_$(problem)_m_$(m)_n_$(n)_idx_$(idx)"
                    sol_path = joinpath(solutions_folder, sol_name)
                    sol_data = matread(sol_path)
                    H = sol_data["H"]
                    H = Matrix(H)
                end

                H_norm_0 = matrix_norm_0(H)
                H_norm_1 = norm(H, 1)
                H_rank = calculate_rank(H)

                H_div_bound_norm_0 = H_norm_0 / (m * r + (m- r) * (n - r))
                H_div_AMP_norm_0 = H_norm_0 / matrix_norm_0(AMP)
                H_div_AMP_norm_1 = H_norm_1 / norm(AMP, 1)
                H_rank_ratio = H_rank / r

                H_norm_0 = matrix_norm_0(H)
                H_norm_1 = norm(H, 1)

                H_div_AMP_norm_0 = H_norm_0 / matrix_norm_0(AMP)
                H_div_AMP_norm_1 = H_norm_1 / norm(AMP, 1)

                push!(H_div_bound_norm_0_list, H_div_bound_norm_0)
                push!(H_div_AMP_norm_0_list, H_div_AMP_norm_0)
                push!(H_div_AMP_norm_1_list, H_div_AMP_norm_1)
                push!(H_rank_ratio_list, H_rank_ratio)
            end
        end
    end

    H_div_bound_norm_0_min = minimum(H_div_bound_norm_0_list)
    H_div_bound_norm_0_mean = mean(H_div_bound_norm_0_list)
    H_div_bound_norm_0_max = maximum(H_div_bound_norm_0_list)

    H_div_AMP_norm_0_min = minimum(H_div_AMP_norm_0_list)
    H_div_AMP_norm_0_mean = mean(H_div_AMP_norm_0_list)
    H_div_AMP_norm_0_max = maximum(H_div_AMP_norm_0_list)

    H_div_AMP_norm_1_min = minimum(H_div_AMP_norm_1_list)
    H_div_AMP_norm_1_mean = mean(H_div_AMP_norm_1_list)
    H_div_AMP_norm_1_max = maximum(H_div_AMP_norm_1_list)

    H_rank_ratio_min = minimum(H_rank_ratio_list)
    H_rank_ratio_mean = mean(H_rank_ratio_list)
    H_rank_ratio_max = maximum(H_rank_ratio_list)

    result = DataFrame(
        method = [method],
        H_div_bound_norm_0_min = [H_div_bound_norm_0_min],
        H_div_bound_norm_0_mean = [H_div_bound_norm_0_mean],
        H_div_bound_norm_0_max = [H_div_bound_norm_0_max],
        H_div_AMP_norm_0_min = [H_div_AMP_norm_0_min],
        H_div_AMP_norm_0_mean = [H_div_AMP_norm_0_mean],
        H_div_AMP_norm_0_max = [H_div_AMP_norm_0_max],
        H_div_AMP_norm_1_min = [H_div_AMP_norm_1_min],
        H_div_AMP_norm_1_mean = [H_div_AMP_norm_1_mean],
        H_div_AMP_norm_1_max = [H_div_AMP_norm_1_max],
        H_rank_ratio_min = [H_rank_ratio_min],
        H_rank_ratio_mean = [H_rank_ratio_mean],
        H_rank_ratio_max = [H_rank_ratio_max],
        m_max = [m_max]
    )

    append!(df, result)
end

results_filename = "results_$(exp)_$(problem)_new_table.csv"
results_filepath = joinpath(results_folder, results_filename)
CSV.write(results_filepath, df)