using CSV
using DataFrames
using MAT
using Base.GC
using Statistics

include("types.jl")
include("utility.jl")
include("solvers.jl")
include("admm.jl")

exp = "7"

matrices_folder = "./Experiment_Matrices/Experiment_" * exp
mat_files = readdir(matrices_folder)

results_folder = "results/Experiment_$exp"

solutions_folder = "./Solutions/Experiment_" * exp

methods = ["Gurobi"]
method = methods[1]

# Gurobi parameters
opt_tol = 10^(-5)
time_limit = 1800

df = DataFrame()

PLS_H_norm_0_list = []
PLS_H_norm_1_list = []
PLS_H_rank_list = []
PLS_time_list = []

PMN_H_norm_0_list = []
PMN_H_norm_1_list = []
PMN_H_rank_list = []
PMN_time_list = []

PMX_H_norm_0_list = []
PMX_H_norm_1_list = []
PMX_H_rank_list = []
PMX_time_list = []

count = 0
min_PLS_unsolvable_m = Inf
min_PMN_unsolvable_m = Inf
min_PMX_unsolvable_m = Inf

for mat_file in mat_files
    global count += 1

    mat_path = joinpath(matrices_folder, mat_file)
    mat_data = matread(mat_path)

    println("Solving for matrix: $mat_path")
    
    A = mat_data["matrix"]
    A = Matrix(A)
    AMP = pinv(A)
    AMP_norm_0 = matrix_norm_0(AMP)
    AMP_norm_1 = norm(AMP, 1)

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

        constraints = ["PLS"]

        GRB_PLS_time = "-"
        GRB_PLS_H_norm_0 = -1.0
        GRB_PLS_H_norm_1 = -1.0
        GRB_PLS_H_rank = -1
        if !("-" in PLS_time_list) && (m < min_PLS_unsolvable_m)
            try
                GRB_PLS_time = @elapsed begin
                    GRB_PLS_H = gurobi_solver(data, constraints, opt_tol, time_limit)
                end
                GRB_PLS_H_norm_0 = matrix_norm_0(GRB_PLS_H)
                GRB_PLS_H_norm_1 = norm(GRB_PLS_H, 1)
                GRB_PLS_H_rank = calculate_rank(GRB_PLS_H)

                solution_filename = "Gurobi/Experiment_$(exp)_PLS_m_$(m)_n_$(n)_idx_$(idx)"
                solution_filepath = joinpath(solutions_folder, solution_filename)
                matwrite(solution_filepath, Dict("H" => GRB_PLS_H, "time" => GRB_PLS_time))
            catch e
                if isa(e, ErrorException)
                    GRB_PLS_time = "-"
                    global min_PLS_unsolvable_m = min(m, min_PLS_unsolvable_m)
                else
                    throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
                end
            end
        end

        push!(PLS_H_norm_0_list, GRB_PLS_H_norm_0 / AMP_norm_0)
        push!(PLS_H_norm_1_list, GRB_PLS_H_norm_1 / AMP_norm_1)
        push!(PLS_H_rank_list, GRB_PLS_H_rank / r)
        push!(PLS_time_list, GRB_PLS_time)

        GC.gc()

        constraints = ["PMN"]

        GRB_PMN_time = "-"
        GRB_PMN_H_norm_0 = -1.0
        GRB_PMN_H_norm_1 = -1.0
        GRB_PMN_H_rank = -1
        if !("-" in PMN_time_list) && (m < min_PMN_unsolvable_m)
            try
                GRB_PMN_time = @elapsed begin
                    GRB_PMN_H = gurobi_solver(data, constraints, opt_tol, time_limit)
                end
                GRB_PMN_H_norm_0 = matrix_norm_0(GRB_PMN_H)
                GRB_PMN_H_norm_1 = norm(GRB_PMN_H, 1)
                GRB_PMN_H_rank = calculate_rank(GRB_PMN_H)

                solution_filename = "Gurobi/Experiment_$(exp)_PMN_m_$(m)_n_$(n)_idx_$(idx)"
                solution_filepath = joinpath(solutions_folder, solution_filename)
                matwrite(solution_filepath, Dict("H" => GRB_PMN_H, "time" => GRB_PMN_time))
            catch e
                if isa(e, ErrorException)
                    GRB_PMN_time = "-"
                    global min_PMN_unsolvable_m = min(m, min_PMN_unsolvable_m)
                else
                    throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
                end
            end
        end

        push!(PMN_H_norm_0_list, GRB_PMN_H_norm_0 / AMP_norm_0)
        push!(PMN_H_norm_1_list, GRB_PMN_H_norm_1 / AMP_norm_1)
        push!(PMN_H_rank_list, GRB_PMN_H_rank / r)
        push!(PMN_time_list, GRB_PMN_time)

        GC.gc()

        constraints = ["PMX"]

        GRB_PMX_time = "-"
        GRB_PMX_H_norm_0 = -1.0
        GRB_PMX_H_norm_1 = -1.0
        GRB_PMX_H_rank = -1
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

        push!(PMX_H_norm_0_list, GRB_PMX_H_norm_0 / AMP_norm_0)
        push!(PMX_H_norm_1_list, GRB_PMX_H_norm_1 / AMP_norm_1)
        push!(PMX_H_rank_list, GRB_PMX_H_rank / r)
        push!(PMX_time_list, GRB_PMX_time)

        GC.gc()

        if count % 3 == 0
            GRB_PLS_H_norm_0_mean = -1.0
            GRB_PLS_H_norm_1_mean = -1.0
            GRB_PLS_H_rank_mean = -1.0
            GRB_PLS_time_mean = -1.0

            if !("-" in PLS_time_list)
                GRB_PLS_H_norm_0_mean = mean(PLS_H_norm_0_list)
                GRB_PLS_H_norm_1_mean = mean(PLS_H_norm_1_list)
                GRB_PLS_H_rank_mean = mean(PLS_H_rank_list)
                GRB_PLS_time_mean = mean(PLS_time_list)
            end

            GRB_PMN_H_norm_0_mean = -1.0
            GRB_PMN_H_norm_1_mean = -1.0
            GRB_PMN_H_rank_mean = -1.0
            GRB_PMN_time_mean = -1.0

            if !("-" in PMN_time_list)
                GRB_PMN_H_norm_0_mean = mean(PMN_H_norm_0_list)
                GRB_PMN_H_norm_1_mean = mean(PMN_H_norm_1_list)
                GRB_PMN_H_rank_mean = mean(PMN_H_rank_list)
                GRB_PMN_time_mean = mean(PMN_time_list)
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
                r = [r],
                GRB_PLS_H_norm_0_mean = [GRB_PLS_H_norm_0_mean],
                GRB_PLS_H_norm_1_mean = [GRB_PLS_H_norm_1_mean],
                GRB_PLS_H_rank_mean = [GRB_PLS_H_rank_mean],
                GRB_PLS_time_mean = [GRB_PLS_time_mean],
                GRB_PMN_H_norm_0_mean = [GRB_PMN_H_norm_0_mean],
                GRB_PMN_H_norm_1_mean = [GRB_PMN_H_norm_1_mean],
                GRB_PMN_H_rank_mean = [GRB_PMN_H_rank_mean],
                GRB_PMN_time_mean = [GRB_PMN_time_mean],
                GRB_PMX_H_norm_0_mean = [GRB_PMX_H_norm_0_mean],
                GRB_PMX_H_norm_1_mean = [GRB_PMX_H_norm_1_mean],
                GRB_PMX_H_rank_mean = [GRB_PMX_H_rank_mean],
                GRB_PMX_time_mean = [GRB_PMX_time_mean]
            )

            append!(df, result)

            empty!(PLS_H_norm_0_list)
            empty!(PLS_H_norm_1_list)
            empty!(PLS_H_rank_list)
            empty!(PLS_time_list)

            empty!(PMN_H_norm_0_list)
            empty!(PMN_H_norm_1_list)
            empty!(PMN_H_rank_list)
            empty!(PMN_time_list)

            empty!(PMX_H_norm_0_list)
            empty!(PMX_H_norm_1_list)
            empty!(PMX_H_rank_list)
            empty!(PMX_time_list)
        end

        GC.gc()
    else
        throw(ErrorException("Invalid method chose."))
    end
end

if method == "Gurobi"
    results_filename = "results_$(exp)_GRB.csv"
    results_filepath = joinpath(results_folder, results_filename)
    CSV.write(results_filepath, df)
else
    throw(ErrorException("Invalid method chose."))
end