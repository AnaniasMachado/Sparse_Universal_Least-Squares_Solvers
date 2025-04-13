using CSV
using DataFrames
using MAT
using Base.GC
using Statistics

include("types.jl")
include("utility.jl")
include("solvers.jl")

exp = "1"

matrices_folder = "./Experiment_Matrices/Experiment_" * exp
mat_files = readdir(matrices_folder)

solutions_folder = "./Solutions/Experiment_" * exp

results_folder = "results/Experiment_$exp"

methods = ["Gurobi"]
method = methods[1]

# Gurobi parameters
opt_tol = 10^(-5)
time_limit = 1200

df = DataFrame()

P1_P3_time_list = []

PLS_H_norm_0_list = []
PLS_H_norm_1_list = []
PLS_H_rank_list = []
PLS_time_list = []

count = 0
min_P1_P3_unsolvable_m = 140
min_PLS_unsolvable_m = Inf

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

    # if m <= 500
    #     continue
    # end

    if method == "Gurobi"
        data = DataInst(A, m, n, r)

        constraints = ["P1", "P3"]

        GRB_P1_P3_time = "-"
        if !("-" in P1_P3_time_list) && (m < min_P1_P3_unsolvable_m)
            try
                GRB_P1_P3_time = @elapsed begin
                    GRB_P1_P3_H = gurobi_solver(data, constraints, opt_tol, time_limit)
                end

                solution_filename = "Gurobi/Experiment_$(exp)_P1_P3_m_$(m)_n_$(n)_idx_$(idx)"
                solution_filepath = joinpath(solutions_folder, solution_filename)
                matwrite(solution_filepath, Dict("H" => GRB_P1_P3_H, "time" => GRB_P1_P3_time))
            catch e
                if isa(e, ErrorException)
                    GRB_P1_P3_time = "-"
                    global min_P1_P3_unsolvable_m = min(m, min_P1_P3_unsolvable_m)
                else
                    throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
                end
            end
        end

        push!(P1_P3_time_list, GRB_P1_P3_time)

        constraints = ["PLS"]

        GRB_PLS_time = "-"
        GRB_PLS_H_norm_0 = "-"
        GRB_PLS_H_norm_1 = "-"
        GRB_PLS_H_rank = "-"
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

        push!(PLS_H_norm_0_list, GRB_PLS_H_norm_0)
        push!(PLS_H_norm_1_list, GRB_PLS_H_norm_1)
        push!(PLS_H_rank_list, GRB_PLS_H_rank)
        push!(PLS_time_list, GRB_PLS_time)

        if count % 5 == 0
            GRB_P1_P3_time_mean = -1.0
            if !("-" in P1_P3_time_list)
                GRB_P1_P3_time_mean = mean(P1_P3_time_list)
            end

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

            result = DataFrame(
                m = [m],
                n = [n],
                r = [r],
                GRB_PLS_H_norm_0_mean = [GRB_PLS_H_norm_0_mean],
                GRB_PLS_H_norm_1_mean = [GRB_PLS_H_norm_1_mean],
                GRB_PLS_H_rank_mean = [GRB_PLS_H_rank_mean],
                GRB_P1_P3_time_mean = [GRB_P1_P3_time_mean],
                GRB_PLS_time_mean = [GRB_PLS_time_mean]
            )

            append!(df, result)

            empty!(P1_P3_time_list)

            empty!(PLS_H_norm_0_list)
            empty!(PLS_H_norm_1_list)
            empty!(PLS_H_rank_list)
            empty!(PLS_time_list)
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