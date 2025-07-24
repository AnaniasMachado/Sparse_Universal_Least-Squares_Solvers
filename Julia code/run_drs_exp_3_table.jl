using CSV
using DataFrames
using MAT
using Statistics

include("utility.jl")

exp = "3"

matrices_folder = "./Experiment_Matrices/DRS_Experiment_$(exp)"
mat_files_exp = "5"

methods = ["ADMM", "ADMMe", "DRS", "DRS_Boyd", "DRS_FP"]
method = methods[5]

problems = ["P134"]
problem = problems[1]

solutions_folder = "./Solutions/DRS_Experiment_$(exp)/$(method)"

results_folder = "results/DRS_Experiment_$exp"

m_values = [60+i*20 for i in 0:27]
d_values = [100]

df = DataFrame()

for m in m_values
    n = m
    r = Int(m / 4)
    for d in d_values
        H_div_mr_norm_0_list = []
        H_div_AMP_norm_0_list = []
        H_div_AMP_norm_1_list = []
        H_rank_ratio_list = []
        time_list = []
        for idx in 1:5
            try
                sol_name = "Experiment_$(exp)_$(problem)_m_$(m)_n_$(n)_idx_5"
                sol_path = joinpath(solutions_folder, sol_name)
                sol_data = matread(sol_path)
            catch
                println("Caught error. m = $(m)")
                push!(H_div_mr_norm_0_list, -1)
                break
            end
            mat_name = "experiment_$(mat_files_exp)_matrix_m$(m)_n$(n)_r$(r)_d$(d)_idx$(idx).mat"
            mat_path = joinpath(matrices_folder, mat_name)
            mat_data = matread(mat_path)
            A = mat_data["matrix"]
            A = Matrix(A)
            AMP = pinv(A)

            sol_name = "Experiment_$(exp)_$(problem)_m_$(m)_n_$(n)_idx_$(idx)"
            sol_path = joinpath(solutions_folder, sol_name)
            sol_data = matread(sol_path)
            H = sol_data["H"]
            H = Matrix(H)
            time = sol_data["time"]

            H_norm_0 = matrix_norm_0(H)
            H_norm_1 = norm(H, 1)
            H_rank = calculate_rank(H)

            H_div_mr_norm_0 = H_norm_0 / (m * r + n * r - r^2)
            H_div_AMP_norm_0 = H_norm_0 / matrix_norm_0(AMP)
            H_div_AMP_norm_1 = H_norm_1 / norm(AMP, 1)
            H_rank_ratio = H_rank / r

            push!(H_div_mr_norm_0_list, H_div_mr_norm_0)
            push!(H_div_AMP_norm_0_list, H_div_AMP_norm_0)
            push!(H_div_AMP_norm_1_list, H_div_AMP_norm_1)
            push!(H_rank_ratio_list, H_rank_ratio)
            push!(time_list, time)
        end

        if !(-1 in H_div_mr_norm_0_list)
            H_div_mr_norm_0_mean = mean(H_div_mr_norm_0_list)
            H_div_AMP_norm_0_mean = mean(H_div_AMP_norm_0_list)
            H_div_AMP_norm_1_mean = mean(H_div_AMP_norm_1_list)
            H_rank_ratio_mean = mean(H_rank_ratio_list)
            time_mean = mean(time_list)

            result = DataFrame(
                m = [m],
                r = [r],
                d = [d],
                H_div_mr_norm_0_mean = [H_div_mr_norm_0_mean],
                H_div_AMP_norm_0_mean = [H_div_AMP_norm_0_mean],
                H_div_AMP_norm_1_mean = [H_div_AMP_norm_1_mean],
                H_rank_ratio_mean = [H_rank_ratio_mean],
                time_mean = [time_mean]
            )

            append!(df, result)
        end
    end
end

results_filename = "results_$(exp)_$(method)_$(problem)_table.csv"
results_filepath = joinpath(results_folder, results_filename)
CSV.write(results_filepath, df)