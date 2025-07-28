using CSV
using DataFrames
using MAT
using Statistics

include("utility.jl")

function count_zeros(H)
    n, m = size(H)
    count = 0
    for j in 1:m
        for i in 1:n
            if H[i, j] == 0
                count += 1
            end
        end
    end
    return count
    return count(==(0), H)
end

function count_between(H)
    count(x -> 10^(-5) < x < 10^(-3), H)
end

exp = "2"
grb_exp = "9"

matrices_folder = "./Experiment_Matrices/DRS_Experiment_$(exp)"
mat_files_exp = "1"

methods = ["Gurobi", "DRS", "DRS_Boyd", "DRS_FP"]

problems = ["P13"]
problem = problems[1]

gurobi_solutions_folder = "./Solutions/Experiment_$(grb_exp)/Gurobi"

results_folder = "results/DRS_Experiment_$exp"

m_values = [60+i*20 for i in 0:27]
d_values = [100]

df = DataFrame()

for method in methods
    solutions_folder = "./Solutions/DRS_Experiment_$(exp)/$(method)"

    H_div_bound_count_0_list = []
    H_div_bound_norm_0_list = []
    H_div_bound_count_between_list = []
    m_max = 600
    for m in m_values
        n = Int(m / 2)
        r = Int(m / 4)
        for d in d_values
            for idx in 1:5
                try
                    if method != "Gurobi"
                        sol_name = "Experiment_$(exp)_$(problem)_m_$(m)_n_$(n)_idx_5"
                        sol_path = joinpath(solutions_folder, sol_name)
                        sol_data = matread(sol_path)
                    else
                        grb_name = "Experiment_$(grb_exp)_PLS_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
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
                    grb_name = "Experiment_$(grb_exp)_PLS_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
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
                H_count_0 = count_zeros(H)
                H_between = count_between(H)

                H_div_bound_norm_0 = H_norm_0 / (m * r)
                H_div_bound_count_0 = H_count_0 / (m * r)
                H_div_bound_count_between = H_between / (m * r)

                push!(H_div_bound_norm_0_list, H_div_bound_norm_0)
                push!(H_div_bound_count_0_list, H_div_bound_count_0)
                push!(H_div_bound_count_between_list, H_div_bound_count_between)
            end
        end
    end

    H_div_bound_norm_0_min = minimum(H_div_bound_norm_0_list)
    H_div_bound_norm_0_mean = mean(H_div_bound_norm_0_list)
    H_div_bound_norm_0_max = maximum(H_div_bound_norm_0_list)

    H_div_bound_count_0_min = minimum(H_div_bound_norm_0_list)
    H_div_bound_count_0_mean = mean(H_div_bound_norm_0_list)
    H_div_bound_count_0_max = maximum(H_div_bound_norm_0_list)

    H_div_bound_count_between_min = minimum(H_div_bound_count_between_list)
    H_div_bound_count_between_mean = mean(H_div_bound_count_between_list)
    H_div_bound_count_between_max = maximum(H_div_bound_count_between_list)

    result = DataFrame(
        method = [method],
        H_div_bound_norm_0_min = [H_div_bound_norm_0_min],
        H_div_bound_norm_0_mean = [H_div_bound_norm_0_mean],
        H_div_bound_norm_0_max = [H_div_bound_norm_0_max],
        H_div_bound_count_0_min = [H_div_bound_count_0_min],
        H_div_bound_count_0_mean = [H_div_bound_count_0_mean],
        H_div_bound_count_0_max = [H_div_bound_count_0_max],
        H_div_bound_count_between_min = [H_div_bound_count_between_min],
        H_div_bound_count_between_mean = [H_div_bound_count_between_mean],
        H_div_bound_count_between_max = [H_div_bound_count_between_max]
    )

    append!(df, result)
end

results_filename = "results_$(exp)_$(problem)_n0_count.csv"
results_filepath = joinpath(results_folder, results_filename)
CSV.write(results_filepath, df)