using CSV
using DataFrames
using MAT
using Base.GC
using Statistics

include("types.jl")
include("utility.jl")
include("solvers.jl")

exp = "3"

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

P123_time_list = []

PLSr_H_norm_0_list = []
PLSr_H_norm_1_list = []
PLSr_H_rank_list = []
PLSr_time_list = []

count = 0
min_P123_unsolvable_m = Inf
min_PLSr_unsolvable_m = Inf

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

    # if m <= 500
    #     continue
    # end

    if method == "Gurobi"
        data = DataInst(A, m, n, r, AMP=AMP)

        constraints = ["P1", "P123", "P3"]

        GRB_P123_time = "-"
        if !("-" in P123_time_list) && (m < min_P123_unsolvable_m)
            try
                GRB_P123_time = @elapsed begin
                    GRB_P123_H = gurobi_solver(data, constraints, opt_tol, time_limit)
                end

                solution_filename = "Gurobi/Experiment_$(exp)_P123_m_$(m)_n_$(n)_idx_$(idx)"
                solution_filepath = joinpath(solutions_folder, solution_filename)
                matwrite(solution_filepath, Dict("H" => GRB_P123_H, "time" => GRB_P123_time))
            catch e
                if isa(e, ErrorException)
                    GRB_P123_time = "-"
                    global min_P123_unsolvable_m = min(m, min_P123_unsolvable_m)
                else
                    throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
                end
            end
        end

        push!(P123_time_list, GRB_P123_time)

        GC.gc()

        constraints = ["PLSr"]

        GRB_PLSr_time = "-"
        GRB_PLSr_H_norm_0 = "-"
        GRB_PLSr_H_norm_1 = "-"
        GRB_PLSr_H_rank = "-"
        if !("-" in PLSr_time_list) && (m < min_PLSr_unsolvable_m)
            try
                GRB_PLSr_time = @elapsed begin
                    GRB_PLSr_H = gurobi_solver(data, constraints, opt_tol, time_limit)
                end
                GRB_PLSr_H_norm_0 = matrix_norm_0(GRB_PLSr_H)
                GRB_PLSr_H_norm_1 = norm(GRB_PLSr_H, 1)
                GRB_PLSr_H_rank = calculate_rank(GRB_PLSr_H)

                solution_filename = "Gurobi/Experiment_$(exp)_PLSr_m_$(m)_n_$(n)_idx_$(idx)"
                solution_filepath = joinpath(solutions_folder, solution_filename)
                matwrite(solution_filepath, Dict("H" => GRB_PLSr_H, "time" => GRB_PLSr_time))
            catch e
                if isa(e, ErrorException)
                    GRB_PLSr_time = "-"
                    global min_PLSr_unsolvable_m = min(m, min_PLSr_unsolvable_m)
                else
                    throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
                end
            end
        end

        push!(PLSr_H_norm_0_list, GRB_PLSr_H_norm_0)
        push!(PLSr_H_norm_1_list, GRB_PLSr_H_norm_1)
        push!(PLSr_H_rank_list, GRB_PLSr_H_rank)
        push!(PLSr_time_list, GRB_PLSr_time)

        if count % 5 == 0
            GRB_P123_time_mean = -1.0
            if !("-" in P123_time_list)
                GRB_P123_time_mean = mean(P123_time_list)
            end

            GRB_PLSr_H_norm_0_mean = -1.0
            GRB_PLSr_H_norm_1_mean = -1.0
            GRB_PLSr_H_rank_mean = -1.0
            GRB_PLSr_time_mean = -1.0

            if !("-" in PLSr_time_list)
                GRB_PLSr_H_norm_0_mean = mean(PLSr_H_norm_0_list)
                GRB_PLSr_H_norm_1_mean = mean(PLSr_H_norm_1_list)
                GRB_PLSr_H_rank_mean = mean(PLSr_H_rank_list)
                GRB_PLSr_time_mean = mean(PLSr_time_list)
            end

            result = DataFrame(
                m = [m],
                n = [n],
                r = [r],
                GRB_PLSr_H_norm_0_mean = [GRB_PLSr_H_norm_0_mean],
                GRB_PLSr_H_norm_1_mean = [GRB_PLSr_H_norm_1_mean],
                GRB_PLSr_H_rank_mean = [GRB_PLSr_H_rank_mean],
                GRB_P123_time_mean = [GRB_P123_time_mean],
                GRB_PLSr_time_mean = [GRB_PLSr_time_mean]
            )

            append!(df, result)

            empty!(P123_time_list)

            empty!(PLSr_H_norm_0_list)
            empty!(PLSr_H_norm_1_list)
            empty!(PLSr_H_rank_list)
            empty!(PLSr_time_list)
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