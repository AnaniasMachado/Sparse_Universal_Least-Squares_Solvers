using CSV
using DataFrames
using MAT
using Base.GC
using Statistics

include("types.jl")
include("utility.jl")
include("solvers.jl")
include("admm.jl")
include("drs.jl")

exp = "1"

matrices_folder = "./Experiment_Matrices/DRS_Experiment_" * exp
mat_files = readdir(matrices_folder)

results_folder = "results/DRS_Experiment_$exp"

solutions_folder = "./Solutions/DRS_Experiment_" * exp

methods = ["ADMM", "DRS"]
method = methods[1]

# Mixed parameters
epsilon = 10^(-5)
eps_abs = epsilon
eps_rel = epsilon
fixed_tol = false
eps_opt = epsilon
time_limit = 1200

# ADMM parameters
rho = 3.0

# DRS paramter
lambda = 10^(-2)
problems = ["P123"]
problem = problems[1]
stop_crits = ["Boyd", "Fixed_Point"]
stop_crit = stop_crits[1]

df = DataFrame()

ADMM_H_norm_0_list = []
ADMM_H_norm_1_list = []
ADMM_H_rank_list = []
ADMM_time_list = []

DRS_H_norm_0_list = []
DRS_H_norm_1_list = []
DRS_H_rank_list = []
DRS_time_list = []

count = 0
min_ADMM_unsolvable_m = Inf
min_DRS_unsolvable_m = Inf

for mat_file in mat_files
    global count += 1

    mat_path = joinpath(matrices_folder, mat_file)
    mat_data = matread(mat_path)

    println("Solving for matrix: $mat_path")
    
    A = mat_data["matrix"]
    A = Matrix(A)

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

    if method == "ADMM"
        ADMM_time = -1.0
        ADMM_H_norm_0 = -1.0
        ADMM_H_norm_1 = -1.0
        ADMM_H_rank = -1.0
        if (m < min_ADMM_unsolvable_m)
            ADMM_time = @elapsed begin
                ADMM_H, ADMM_k = admm_p123(A, rho, eps_abs, eps_rel, fixed_tol, eps_opt, time_limit)
            end
            if ADMM_H == "-"
                global min_ADMM_unsolvable_m = min(m, min_ADMM_unsolvable_m)
            else
                ADMM_H_norm_0 = matrix_norm_0(ADMM_H)
                ADMM_H_norm_1 = norm(ADMM_H, 1)
                ADMM_H_rank = calculate_rank(ADMM_H)

                if fixed_tol
                    solution_filename = "ADMMe/Experiment_$(exp)_P123_m_$(m)_n_$(n)_idx_$(idx)"
                    solution_filepath = joinpath(solutions_folder, solution_filename)
                    matwrite(solution_filepath, Dict("H" => ADMM_H, "time" => ADMM_time, "k" => ADMM_k))
                else
                    solution_filename = "ADMM/Experiment_$(exp)_P123_m_$(m)_n_$(n)_idx_$(idx)"
                    solution_filepath = joinpath(solutions_folder, solution_filename)
                    matwrite(solution_filepath, Dict("H" => ADMM_H, "time" => ADMM_time, "k" => ADMM_k))
                end
            end
        end

        push!(ADMM_H_norm_0_list, ADMM_H_norm_0)
        push!(ADMM_H_norm_1_list, ADMM_H_norm_1)
        push!(ADMM_H_rank_list, ADMM_H_rank)
        push!(ADMM_time_list, ADMM_time)
    
        GC.gc()

        if count % 5 == 0
            ADMM_H_norm_0_mean = -1.0
            ADMM_H_norm_1_mean = -1.0
            ADMM_H_rank_mean = -1.0
            ADMM_time_mean = -1.0

            if !("-" in ADMM_time_list)
                ADMM_H_norm_0_mean = mean(ADMM_H_norm_0_list)
                ADMM_H_norm_1_mean = mean(ADMM_H_norm_1_list)
                ADMM_H_rank_mean = mean(ADMM_H_rank_list)
                ADMM_time_mean = mean(ADMM_time_list)
            end

            result = DataFrame(
                m = [m],
                n = [n],
                r = [r],
                ADMM_H_norm_0_mean = [ADMM_H_norm_0_mean],
                ADMM_H_norm_1_mean = [ADMM_H_norm_1_mean],
                ADMM_H_rank_mean = [ADMM_H_rank_mean],
                ADMM_time_mean = [ADMM_time_mean]
            )

            append!(df, result)

            empty!(ADMM_H_norm_0_list)
            empty!(ADMM_H_norm_1_list)
            empty!(ADMM_H_rank_list)
            empty!(ADMM_time_list)

            GC.gc()
        end
    elseif method == "DRS"
        DRS_time = -1.0
        DRS_H_norm_0 = -1.0
        DRS_H_norm_1 = -1.0
        DRS_H_rank = -1.0
        if (m < min_DRS_unsolvable_m)
            DRS_time = @elapsed begin
                DRS_H, DRS_k = drs(A, lambda, eps_abs, eps_rel, problem, fixed_tol, eps_opt, stop_crit, time_limit)
            end
            if DRS_H == "-"
                global min_DRS_unsolvable_m = min(m, min_DRS_unsolvable_m)
            else
                DRS_H_norm_0 = matrix_norm_0(DRS_H)
                DRS_H_norm_1 = norm(DRS_H, 1)
                DRS_H_rank = calculate_rank(DRS_H)

                if fixed_tol
                    solution_filename = "DRS/Experiment_$(exp)_P123_m_$(m)_n_$(n)_idx_$(idx)"
                    solution_filepath = joinpath(solutions_folder, solution_filename)
                    matwrite(solution_filepath, Dict("H" => DRS_H, "time" => DRS_time, "k" => DRS_k))
                elseif stop_crit == "Boyd"
                    solution_filename = "DRS_Boyd/Experiment_$(exp)_P123_m_$(m)_n_$(n)_idx_$(idx)"
                    solution_filepath = joinpath(solutions_folder, solution_filename)
                    matwrite(solution_filepath, Dict("H" => DRS_H, "time" => DRS_time, "k" => DRS_k))
                elseif stop_crit == "Fixed_Point"
                    solution_filename = "DRS_FP/Experiment_$(exp)_P123_m_$(m)_n_$(n)_idx_$(idx)"
                    solution_filepath = joinpath(solutions_folder, solution_filename)
                    matwrite(solution_filepath, Dict("H" => DRS_H, "time" => DRS_time, "k" => DRS_k))
                end
            end
        end

        push!(DRS_H_norm_0_list, DRS_H_norm_0)
        push!(DRS_H_norm_1_list, DRS_H_norm_1)
        push!(DRS_H_rank_list, DRS_H_rank)
        push!(DRS_time_list, DRS_time)
    
        GC.gc()

        if count % 5 == 0
            DRS_H_norm_0_mean = -1.0
            DRS_H_norm_1_mean = -1.0
            DRS_H_rank_mean = -1.0
            DRS_time_mean = -1.0

            if !("-" in DRS_time_list)
                DRS_H_norm_0_mean = mean(DRS_H_norm_0_list)
                DRS_H_norm_1_mean = mean(DRS_H_norm_1_list)
                DRS_H_rank_mean = mean(DRS_H_rank_list)
                DRS_time_mean = mean(DRS_time_list)
            end

            result = DataFrame(
                m = [m],
                n = [n],
                r = [r],
                DRS_H_norm_0_mean = [DRS_H_norm_0_mean],
                DRS_H_norm_1_mean = [DRS_H_norm_1_mean],
                DRS_H_rank_mean = [DRS_H_rank_mean],
                DRS_time_mean = [DRS_time_mean]
            )

            append!(df, result)

            empty!(DRS_H_norm_0_list)
            empty!(DRS_H_norm_1_list)
            empty!(DRS_H_rank_list)
            empty!(DRS_time_list)

            GC.gc()
        end
    else
        throw(ErrorException("Invalid method chose."))
    end
end

if method == "ADMM"
    if fixed_tol
        results_filename = "results_$(exp)_ADMMe.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    else
        results_filename = "results_$(exp)_ADMM.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    end
elseif method == "DRS"
    if fixed_tol
        results_filename = "results_$(exp)_DRS.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    elseif !fixed_tol && stop_crit == "Boyd"
        results_filename = "results_$(exp)_DRS_Boyd.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    else
        results_filename = "results_$(exp)_DRS_FP.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    end
else
    throw(ErrorException("Invalid method chose."))
end