using CSV
using DataFrames
using MAT
using Base.GC
using Statistics

include("types.jl")
include("utility.jl")
include("solvers.jl")
include("drs.jl")

exp = "8"

matrices_folder = "./Experiment_Matrices/Experiment_" * exp
mat_files = readdir(matrices_folder)

results_folder = "results/Experiment_$exp"

solutions_folder = "./Solutions/Experiment_" * exp

methods = ["DRS"]
method = methods[1]

# ADMM parameters
rho = 3.0

# DRS paramter
lambda = 10^(-2)
problems = ["PLS", "PMN", "P134", "P123", "P124"]
problem = problems[5]
stop_crits = ["Fixed_Point"]
stop_crit = stop_crits[1]

epsilon = 10^(-5)
eps_abs = epsilon
eps_rel = epsilon
fixed_tol = false
eps_opt = epsilon
time_limit = 1800

df = DataFrame()

DRS_H_div_AMP_norm_0_list = []
DRS_H_div_AMP_norm_1_list = []
DRS_H_rank_ratio_list = []
DRS_time_list = []

count = 0
num_instances = 3
min_DRS_unsolvable_m = Dict()

for mat_file in mat_files
    global count += 1

    mat_path = joinpath(matrices_folder, mat_file)
    mat_data = matread(mat_path)

    println("Solving for matrix: $mat_path")
    
    A = mat_data["matrix"]
    A = Matrix(A)
    AMP = pinv(A)

    m_value = match(r"m(\d+)", mat_file).captures[1]
    n_value = match(r"n(\d+)", mat_file).captures[1]
    r_value = match(r"r(\d+)", mat_file).captures[1]
    d_value = match(r"d(\d+)", mat_file).captures[1]
    idx_value = match(r"idx(\d+)", mat_file).captures[1]

    m = parse(Int, m_value)
    n = parse(Int, n_value)
    r = parse(Int, r_value)
    d = parse(Int, d_value)
    idx = parse(Int, idx_value)

    if !haskey(min_DRS_unsolvable_m, d)
        min_DRS_unsolvable_m[d] = Inf
    end

    if problem == "P124"
        t = m
        m = n
        n = t
    end

    if method == "DRS"
        DRS_time = -1.0
        DRS_H_norm_0 = -1.0
        DRS_H_norm_1 = -1.0
        DRS_H_rank = -1.0
        if (m < min_DRS_unsolvable_m[d])
            DRS_time = @elapsed begin
                DRS_H, DRS_k = drs(A, lambda, eps_abs, eps_rel, problem, fixed_tol, eps_opt, stop_crit, time_limit)
            end
            if DRS_H == "-"
                global min_DRS_unsolvable_m[d] = min(m, min_DRS_unsolvable_m[d])
            else
                DRS_H_norm_0 = matrix_norm_0(DRS_H)
                DRS_H_norm_1 = norm(DRS_H, 1)
                DRS_H_rank = calculate_rank(DRS_H)

                DRS_H_norm_0 = DRS_H_norm_0 / matrix_norm_0(AMP)
                DRS_H_norm_1 = DRS_H_norm_1 / norm(AMP, 1)
                DRS_H_rank = (DRS_H_rank - r) / (m - r)

                if fixed_tol
                    solution_filename = "DRS/Experiment_$(exp)_$(problem)_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                    solution_filepath = joinpath(solutions_folder, solution_filename)
                    matwrite(solution_filepath, Dict("H" => DRS_H, "time" => DRS_time, "k" => DRS_k))
                elseif stop_crit == "Boyd"
                    solution_filename = "DRS_Boyd/Experiment_$(exp)_$(problem)_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                    solution_filepath = joinpath(solutions_folder, solution_filename)
                    matwrite(solution_filepath, Dict("H" => DRS_H, "time" => DRS_time, "k" => DRS_k))
                elseif stop_crit == "Fixed_Point"
                    solution_filename = "DRS_FP/Experiment_$(exp)_$(problem)_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                    solution_filepath = joinpath(solutions_folder, solution_filename)
                    matwrite(solution_filepath, Dict("H" => DRS_H, "time" => DRS_time, "k" => DRS_k))
                end
            end
        end

        push!(DRS_H_div_AMP_norm_0_list, DRS_H_norm_0)
        push!(DRS_H_div_AMP_norm_1_list, DRS_H_norm_1)
        push!(DRS_H_rank_ratio_list, DRS_H_rank)
        push!(DRS_time_list, DRS_time)
    
        GC.gc()

        if count % num_instances == 0
            DRS_H_div_AMP_norm_0_mean = -1.0
            DRS_H_div_AMP_norm_1_mean = -1.0
            DRS_H_rank_ratio_mean = -1.0
            DRS_time_mean = -1.0

            if !(-1.0 in DRS_H_div_AMP_norm_0_list)
                DRS_H_div_AMP_norm_0_mean = mean(DRS_H_div_AMP_norm_0_list)
                DRS_H_div_AMP_norm_1_mean = mean(DRS_H_div_AMP_norm_1_list)
                DRS_H_rank_ratio_mean = mean(DRS_H_rank_ratio_list)
                DRS_time_mean = mean(DRS_time_list)
            end

            result = DataFrame(
                m = [m],
                r = [r],
                d = [d],
                A_norm_0 = [matrix_norm_0(A)],
                A_norm_1 = [norm(A, 1)],
                AMP_norm_0 = [matrix_norm_0(AMP)],
                AMP_norm_1 = [norm(AMP, 1)],
                DRS_H_div_AMP_norm_0_mean = [DRS_H_div_AMP_norm_0_mean],
                DRS_H_div_AMP_norm_1_mean = [DRS_H_div_AMP_norm_1_mean],
                DRS_H_rank_ratio_mean = [DRS_H_rank_ratio_mean],
                DRS_time_mean = [DRS_time_mean]
            )

            append!(df, result)

            empty!(DRS_H_div_AMP_norm_0_list)
            empty!(DRS_H_div_AMP_norm_1_list)
            empty!(DRS_H_rank_ratio_list)
            empty!(DRS_time_list)

            GC.gc()
        end
    else
        throw(ErrorException("Invalid method chose."))
    end
end

if method == "DRS"
    if fixed_tol
        results_filename = "results_$(exp)_DRS_$(problem).csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    elseif !fixed_tol && stop_crit == "Boyd"
        results_filename = "results_$(exp)_DRS_Boyd_$(problem).csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    else
        results_filename = "results_$(exp)_DRS_FP_$(problem).csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    end
else
    throw(ErrorException("Invalid method chose."))
end