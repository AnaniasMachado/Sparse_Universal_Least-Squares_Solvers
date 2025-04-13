using CSV
using DataFrames
using MAT
using Base.GC
using Statistics

include("types.jl")
include("utility.jl")
include("solvers.jl")

exp = "4"

matrices_folder = "./Experiment_Matrices/Experiment_" * exp
mat_files = readdir(matrices_folder)

results_folder = "results/Experiment_$exp"

solutions_folder = "./Solutions/Experiment_" * exp

methods = ["Gurobi"]
method = methods[1]

# Gurobi parameters
opt_tol = 10^(-5)
time_limit = 1200

df = DataFrame()

P124_time_list = []

PMNr_H_norm_0_list = []
PMNr_H_norm_1_list = []
PMNr_H_rank_list = []
PMNr_time_list = []

count = 0
min_P124_unsolvable_m = Inf
min_PMNr_unsolvable_m = Inf

for mat_file in mat_files
    global count += 1

    mat_path = joinpath(matrices_folder, mat_file)
    mat_data = matread(mat_path)

    println("Solving for matrix: $mat_path")
    
    A = mat_data["matrix"]
    A = Matrix(A)
    A = Matrix(A')
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

    # if m <= 500
    #     continue
    # end

    if method == "Gurobi"
        data = DataInst(A, n, m, r, AMP=AMP)

        constraints = ["P1", "P124", "P4"]

        GRB_P124_time = "-"
        if !("-" in P124_time_list) && (m < min_P124_unsolvable_m)
            try
                GRB_P124_time = @elapsed begin
                    GRB_P124_H = gurobi_solver(data, constraints, opt_tol, time_limit)
                end

                solution_filename = "Gurobi/Experiment_$(exp)_P124_m_$(m)_n_$(n)_idx_$(idx)"
                solution_filepath = joinpath(solutions_folder, solution_filename)
                matwrite(solution_filepath, Dict("H" => GRB_P124_H, "time" => GRB_P124_time))
            catch e
                if isa(e, ErrorException)
                    GRB_P124_time = "-"
                    global min_P124_unsolvable_m = min(m, min_P124_unsolvable_m)
                else
                    throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
                end
            end
        end

        push!(P124_time_list, GRB_P124_time)

        GC.gc()

        constraints = ["PMNr"]

        GRB_PMNr_time = "-"
        GRB_PMNr_H_norm_0 = "-"
        GRB_PMNr_H_norm_1 = "-"
        GRB_PMNr_H_rank = "-"
        if !("-" in PMNr_time_list) && (m < min_PMNr_unsolvable_m)
            try
                GRB_PMNr_time = @elapsed begin
                    GRB_PMNr_H = gurobi_solver(data, constraints, opt_tol, time_limit)
                end
                GRB_PMNr_H_norm_0 = matrix_norm_0(GRB_PMNr_H)
                GRB_PMNr_H_norm_1 = norm(GRB_PMNr_H, 1)
                GRB_PMNr_H_rank = calculate_rank(GRB_PMNr_H)

                solution_filename = "Gurobi/Experiment_$(exp)_PMNr_m_$(m)_n_$(n)_idx_$(idx)"
                solution_filepath = joinpath(solutions_folder, solution_filename)
                matwrite(solution_filepath, Dict("H" => GRB_PMNr_H, "time" => GRB_PMNr_time))
            catch e
                if isa(e, ErrorException)
                    GRB_PMNr_time = "-"
                    global min_PMNr_unsolvable_m = min(m, min_PMNr_unsolvable_m)
                else
                    throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
                end
            end
        end

        push!(PMNr_H_norm_0_list, GRB_PMNr_H_norm_0)
        push!(PMNr_H_norm_1_list, GRB_PMNr_H_norm_1)
        push!(PMNr_H_rank_list, GRB_PMNr_H_rank)
        push!(PMNr_time_list, GRB_PMNr_time)

        if count % 5 == 0
            GRB_P124_time_mean = -1.0
            if !("-" in P124_time_list)
                GRB_P124_time_mean = mean(P124_time_list)
            end

            GRB_PMNr_H_norm_0_mean = -1.0
            GRB_PMNr_H_norm_1_mean = -1.0
            GRB_PMNr_H_rank_mean = -1.0
            GRB_PMNr_time_mean = -1.0

            if !("-" in PMNr_time_list)
                GRB_PMNr_H_norm_0_mean = mean(PMNr_H_norm_0_list)
                GRB_PMNr_H_norm_1_mean = mean(PMNr_H_norm_1_list)
                GRB_PMNr_H_rank_mean = mean(PMNr_H_rank_list)
                GRB_PMNr_time_mean = mean(PMNr_time_list)
            end

            result = DataFrame(
                m = [n],
                n = [m],
                r = [r],
                GRB_PLSr_H_norm_0_mean = [GRB_PMNr_H_norm_0_mean],
                GRB_PLSr_H_norm_1_mean = [GRB_PMNr_H_norm_1_mean],
                GRB_PLSr_H_rank_mean = [GRB_PMNr_H_rank_mean],
                GRB_P123_time_mean = [GRB_P124_time_mean],
                GRB_PLSr_time_mean = [GRB_PMNr_time_mean]
            )

            append!(df, result)

            empty!(P124_time_list)

            empty!(PMNr_H_norm_0_list)
            empty!(PMNr_H_norm_1_list)
            empty!(PMNr_H_rank_list)
            empty!(PMNr_time_list)
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