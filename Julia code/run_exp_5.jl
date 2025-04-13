using CSV
using DataFrames
using MAT
using Base.GC
using Statistics

include("types.jl")
include("utility.jl")
include("solvers.jl")
include("admm.jl")

exp = "5"

matrices_folder = "./Experiment_Matrices/Experiment_" * exp
mat_files = readdir(matrices_folder)

results_folder = "results/Experiment_$exp"

solutions_folder = "./Solutions/Experiment_" * exp

methods = ["Gurobi", "ADMM"]
method = methods[2]

# Gurobi parameters
opt_tol = 10^(-5)

# ADMM parameters
rho = 3.0
epsilon = 10^(-5)
eps_abs = epsilon
eps_rel = epsilon
fixed_tol = true
eps_opt = epsilon
time_limit = 1200

df = DataFrame()

P1_P3_P4_time_list = []
PLS_P4_time_list = []
PMN_P3_time_list = []

PMX_H_norm_0_list = []
PMX_H_norm_1_list = []
PMX_H_rank_list = []
PMX_time_list = []

ADMM_H_norm_0_list = []
ADMM_H_norm_1_list = []
ADMM_H_rank_list = []
ADMM_time_list = []

ADMMe_H_norm_0_list = []
ADMMe_H_norm_1_list = []
ADMMe_H_rank_list = []
ADMMe_time_list = []

count = 0
min_P1_P3_P4_unsolvable_m = Inf
min_PLS_P4_unsolvable_m = Inf
min_PMN_P3_unsolvable_m = Inf
min_PMX_unsolvable_m = Inf
min_ADMM_unsolvable_m = Inf
min_ADMMe_unsolvable_m = 100

for mat_file in mat_files
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

    if method == "Gurobi"
        data = DataInst(A, m, n, r)

        constraints = ["P1", "P3", "P4"]

        GRB_P1_P3_P4_time = "-"
        if !("-" in P1_P3_P4_time_list) && (m < min_P1_P3_P4_unsolvable_m)
            try
                GRB_P1_P3_P4_time = @elapsed begin
                    GRB_P1_P3_P4_H = gurobi_solver(data, constraints, opt_tol, time_limit)
                end

                solution_filename = "Gurobi/Experiment_$(exp)_P1_P3_P4_m_$(m)_n_$(n)_idx_$(idx)"
                solution_filepath = joinpath(solutions_folder, solution_filename)
                matwrite(solution_filepath, Dict("H" => GRB_P1_P3_P4_H, "time" => GRB_P1_P3_P4_time))
            catch e
                if isa(e, ErrorException)
                    GRB_P1_P3_P4_time = "-"
                    global min_P1_P3_P4_unsolvable_m = min(m, min_P1_P3_P4_unsolvable_m)
                else
                    throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
                end
            end
        end

        push!(P1_P3_P4_time_list, GRB_P1_P3_P4_time)

        GC.gc()

        constraints = ["PLS", "P4"]

        GRB_PLS_P4_time = "-"
        if !("-" in PLS_P4_time_list) && (m < min_PLS_P4_unsolvable_m)
            try
                GRB_PLS_P4_time = @elapsed begin
                    GRB_PLS_P4_H = gurobi_solver(data, constraints, opt_tol, time_limit)
                end

                solution_filename = "Gurobi/Experiment_$(exp)_PLS_P4_m_$(m)_n_$(n)_idx_$(idx)"
                solution_filepath = joinpath(solutions_folder, solution_filename)
                matwrite(solution_filepath, Dict("H" => GRB_PLS_P4_H, "time" => GRB_PLS_P4_time))
            catch e
                if isa(e, ErrorException)
                    GRB_PLS_P4_time = "-"
                    global min_PLS_P4_unsolvable_m = min(m, min_PLS_P4_unsolvable_m)
                else
                    throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
                end
            end
        end

        push!(PLS_P4_time_list, GRB_PLS_P4_time)

        GC.gc()

        constraints = ["PMN", "P3"]

        GRB_PMN_P3_time = "-"
        if !("-" in PMN_P3_time_list) && (m < min_PMN_P3_unsolvable_m)
            try
                GRB_PMN_P3_time = @elapsed begin
                    GRB_PMN_P3_H = gurobi_solver(data, constraints, opt_tol, time_limit)
                end

                solution_filename = "Gurobi/Experiment_$(exp)_PMN_P3_m_$(m)_n_$(n)_idx_$(idx)"
                solution_filepath = joinpath(solutions_folder, solution_filename)
                matwrite(solution_filepath, Dict("H" => GRB_PMN_P3_H, "time" => GRB_PMN_P3_time))
            catch e
                if isa(e, ErrorException)
                    GRB_PMN_P3_time = "-"
                    global min_PMN_P3_unsolvable_m = min(m, min_PMN_P3_unsolvable_m)
                else
                    throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
                end
            end
        end

        push!(PMN_P3_time_list, GRB_PMN_P3_time)

        GC.gc()

        constraints = ["PMX"]

        GRB_PMX_time = "-"
        GRB_PMX_H_norm_0 = "-"
        GRB_PMX_H_norm_1 = "-"
        GRB_PMX_H_rank = "-"
        if !("-" in PMX_time_list) && (m < min_PMX_unsolvable_m)
            try
                GRB_PMX_time = @elapsed begin
                    GRB_PMX_H = gurobi_solver(data, constraints, opt_tol, time_limit)
                end
                GRB_PMX_H_norm_0 = matrix_norm_0(GRB_PMX_H)
                GRB_PMX_H_norm_1 = norm(GRB_PMX_H, 1)
                GRB_PMX_H_rank = calculate_rank(GRB_PMX_H)

                solution_filename = "Gurobi/Experiment_$(exp)_PMX_m_$(m)_n_$(n)_idx_$(idx)"
                solution_filepath = joinpath(solutions_folder, solution_filename)
                matwrite(solution_filepath, Dict("H" => GRB_PMX_H, "time" => GRB_PMX_time))
            catch e
                if isa(e, ErrorException)
                    GRB_PMX_time = "-"
                    global min_PMX_unsolvable_m = min(m, min_PMX_unsolvable_m)
                else
                    throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
                end
            end
        end

        push!(PMX_H_norm_0_list, GRB_PMX_H_norm_0)
        push!(PMX_H_norm_1_list, GRB_PMX_H_norm_1)
        push!(PMX_H_rank_list, GRB_PMX_H_rank)
        push!(PMX_time_list, GRB_PMX_time)
    
        GC.gc()

        if count % 5 == 0
            GRB_P1_P3_P4_time_mean = -1.0
            if !("-" in P1_P3_P4_time_list)
                GRB_P1_P3_P4_time_mean = mean(P1_P3_P4_time_list)
            end

            GRB_PLS_P4_time_mean = -1.0
            if !("-" in PLS_P4_time_list)
                GRB_PLS_P4_time_mean = mean(PLS_P4_time_list)
            end

            GRB_PMN_P3_time_mean = -1.0
            if !("-" in PMN_P3_time_list)
                GRB_PMN_P3_time_mean = mean(PMN_P3_time_list)
            end

            GRB_PMX_H_norm_0_mean = -1.0
            GRB_PMX_H_norm_1_mean = -1.0
            GRB_PMX_H_rank_mean = -1.0
            GRB_PMX_time_mean = -1.0

            if !("-" in PMX_time_list)
                GRB_PMX_H_norm_0_mean = mean(PMX_H_norm_0_list)
                GRB_PMX_H_norm_1_mean = mean(PMX_H_norm_1_list)
                GRB_PMX_H_rank_mean = mean(PMX_H_rank_list)
                GRB_PMX_time_mean = mean(PMX_time_list)
            end

            result = DataFrame(
                m = [m],
                n = [n],
                r = [r],
                GRB_PMX_H_norm_0_mean = [GRB_PMX_H_norm_0_mean],
                GRB_PMX_H_norm_1_mean = [GRB_PMX_H_norm_1_mean],
                GRB_PMX_H_rank_mean = [GRB_PMX_H_rank_mean],
                GRB_P1_P3_P4_time_mean = [GRB_P1_P3_P4_time_mean],
                GRB_PLS_P4_time_mean = [GRB_PLS_P4_time_mean],
                GRB_PMN_P3_time_mean = [GRB_PMN_P3_time_mean],
                GRB_PMX_time_mean = [GRB_PMX_time_mean]
            )

            append!(df, result)

            empty!(P1_P3_P4_time_list)
            empty!(PLS_P4_time_list)
            empty!(PMN_P3_time_list)

            empty!(PMX_H_norm_0_list)
            empty!(PMX_H_norm_1_list)
            empty!(PMX_H_rank_list)
            empty!(PMX_time_list)
        end

        GC.gc()
    elseif method == "ADMM"
        ADMM_time = -1.0
        ADMM_H_norm_0 = -1.0
        ADMM_H_norm_1 = -1.0
        ADMM_H_rank = -1.0
        if (m < min(min_ADMM_unsolvable_m, min_ADMMe_unsolvable_m))
            ADMM_time = @elapsed begin
                ADMM_H = admm_p134(A, rho, eps_abs, eps_rel, fixed_tol, eps_opt, time_limit)
            end
            if fixed_tol
                if ADMM_H == "-"
                    global min_ADMMe_unsolvable_m = min(m, min_ADMMe_unsolvable_m)
                else
                    ADMM_H_norm_0 = matrix_norm_0(ADMM_H)
                    ADMM_H_norm_1 = norm(ADMM_H, 1)
                    ADMM_H_rank = calculate_rank(ADMM_H)

                    solution_filename = "ADMMe/Experiment_$(exp)_P134_m_$(m)_n_$(n)_idx_$(idx)"
                    solution_filepath = joinpath(solutions_folder, solution_filename)
                    matwrite(solution_filepath, Dict("H" => ADMM_H, "time" => ADMM_time))
                end
            else
                if ADMM_H == "-"
                    global min_ADMM_unsolvable_m = min(m, min_ADMM_unsolvable_m)
                else
                    ADMM_H_norm_0 = matrix_norm_0(ADMM_H)
                    ADMM_H_norm_1 = norm(ADMM_H, 1)
                    ADMM_H_rank = calculate_rank(ADMM_H)

                    solution_filename = "ADMM/Experiment_$(exp)_P134_m_$(m)_n_$(n)_idx_$(idx)"
                    solution_filepath = joinpath(solutions_folder, solution_filename)
                    matwrite(solution_filepath, Dict("H" => ADMM_H, "time" => ADMM_time))
                end
            end
        end

        if fixed_tol
            push!(ADMMe_H_norm_0_list, ADMM_H_norm_0)
            push!(ADMMe_H_norm_1_list, ADMM_H_norm_1)
            push!(ADMMe_H_rank_list, ADMM_H_rank)
            push!(ADMMe_time_list, ADMM_time)
        else
            push!(ADMM_H_norm_0_list, ADMM_H_norm_0)
            push!(ADMM_H_norm_1_list, ADMM_H_norm_1)
            push!(ADMM_H_rank_list, ADMM_H_rank)
            push!(ADMM_time_list, ADMM_time)
        end
    
        GC.gc()

        if count % 5 == 0
            ADMM_H_norm_0_mean = -1.0
            ADMM_H_norm_1_mean = -1.0
            ADMM_H_rank_mean = -1.0
            ADMM_time_mean = -1.0

            if fixed_tol
                if !("-" in ADMMe_time_list)
                    ADMM_H_norm_0_mean = mean(ADMMe_H_norm_0_list)
                    ADMM_H_norm_1_mean = mean(ADMMe_H_norm_1_list)
                    ADMM_H_rank_mean = mean(ADMMe_H_rank_list)
                    ADMM_time_mean = mean(ADMMe_time_list)
                end
            else
                if !("-" in ADMM_time_list)
                    ADMM_H_norm_0_mean = mean(ADMM_H_norm_0_list)
                    ADMM_H_norm_1_mean = mean(ADMM_H_norm_1_list)
                    ADMM_H_rank_mean = mean(ADMM_H_rank_list)
                    ADMM_time_mean = mean(ADMM_time_list)
                end
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

            if fixed_tol
                empty!(ADMMe_H_norm_0_list)
                empty!(ADMMe_H_norm_1_list)
                empty!(ADMMe_H_rank_list)
                empty!(ADMMe_time_list)
            else
                empty!(ADMM_H_norm_0_list)
                empty!(ADMM_H_norm_1_list)
                empty!(ADMM_H_rank_list)
                empty!(ADMM_time_list)
            end

            GC.gc()
        end
    else
        throw(ErrorException("Invalid method chose."))
    end
end

# if method == "Gurobi"
#     results_filename = "results_$(exp)_GRB.csv"
#     results_filepath = joinpath(results_folder, results_filename)
#     CSV.write(results_filepath, df)
# elseif method == "ADMM"
#     if fixed_tol
#         results_filename = "results_$(exp)_ADMMe.csv"
#         results_filepath = joinpath(results_folder, results_filename)
#         CSV.write(results_filepath, df)
#     else
#         results_filename = "results_$(exp)_ADMM.csv"
#         results_filepath = joinpath(results_folder, results_filename)
#         CSV.write(results_filepath, df)
#     end
# else
#     throw(ErrorException("Invalid method chose."))
# end