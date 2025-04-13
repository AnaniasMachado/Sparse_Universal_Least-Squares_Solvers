using CSV
using DataFrames
using MAT
using Base.GC
using Statistics

include("types.jl")
include("utility.jl")
include("solvers.jl")

exp = "2"

matrices_folder = "./Experiment_Matrices/Experiment_" * exp
mat_files = readdir(matrices_folder)

results_folder = "results/Experiment_$exp"

solutions_folder = "./Solutions/Experiment_" * exp

results_folder = "results/Experiment_$exp"

methods = ["Gurobi"]
method = methods[1]

# Gurobi parameters
opt_tol = 10^(-5)
time_limit = 1200

df = DataFrame()

P1_P4_time_list = []

PMN_H_norm_0_list = []
PMN_H_norm_1_list = []
PMN_H_rank_list = []
PMN_time_list = []

count = 0
min_P1_P4_unsolvable_m = 60
min_PMN_unsolvable_m = Inf

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

    # if m <= 120
    #     continue
    # end

    if method == "Gurobi"
        data = DataInst(A', n, m, r)

        constraints = ["P1", "P4"]

        GRB_P1_P4_time = "-"
        if !("-" in P1_P4_time_list) && (m < min_P1_P4_unsolvable_m)
            try
                GRB_P1_P4_time = @elapsed begin
                    GRB_P1_P4_H = gurobi_solver(data, constraints, opt_tol, time_limit)
                end

                solution_filename = "Gurobi/Experiment_$(exp)_P1_P4_m_$(m)_n_$(n)_idx_$(idx))"
                solution_filepath = joinpath(solutions_folder, solution_filename)
                matwrite(solution_filepath, Dict("H" => GRB_P1_P4_H, "time" => GRB_P1_P4_time))
            catch e
                if isa(e, ErrorException)
                    GRB_P1_P4_time = "-"
                    global min_P1_P4_unsolvable_m = min(m, min_P1_P4_unsolvable_m)
                else
                    throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
                end
            end
        end

        push!(P1_P4_time_list, GRB_P1_P4_time)

        constraints = ["PMN"]

        GRB_PMN_time = "-"
        GRB_PMN_H_norm_0 = "-"
        GRB_PMN_H_norm_1 = "-"
        GRB_PMN_H_rank = "-"
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

        # GRB_PMN_time = @elapsed begin
        #     GRB_PMN_H = gurobi_solver(data, constraints, opt_tol)
        # end
        # GRB_PMN_H_norm_0 = matrix_norm_0(GRB_PMN_H)
        # GRB_PMN_H_norm_1 = norm(GRB_PMN_H, 1)
        # GRB_PMN_H_rank = calculate_rank(GRB_PMN_H)

        push!(PMN_H_norm_0_list, GRB_PMN_H_norm_0)
        push!(PMN_H_norm_1_list, GRB_PMN_H_norm_1)
        push!(PMN_H_rank_list, GRB_PMN_H_rank)
        push!(PMN_time_list, GRB_PMN_time)

        if count % 5 == 0
            GRB_P1_P4_time_mean = -1.0
            if !("-" in P1_P4_time_list)
                GRB_P1_P4_time_mean = mean(P1_P4_time_list)
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

            result = DataFrame(
                m = [n],
                n = [m],
                r = [r],
                GRB_PMN_H_norm_0_mean = [GRB_PMN_H_norm_0_mean],
                GRB_PMN_H_norm_1_mean = [GRB_PMN_H_norm_1_mean],
                GRB_PMN_H_rank_mean = [GRB_PMN_H_rank_mean],
                GRB_P1_P4_time_mean = [GRB_P1_P4_time_mean],
                GRB_PMN_time_mean = [GRB_PMN_time_mean]
            )

            append!(df, result)

            empty!(P1_P4_time_list)

            empty!(PMN_H_norm_0_list)
            empty!(PMN_H_norm_1_list)
            empty!(PMN_H_rank_list)
            empty!(PMN_time_list)
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