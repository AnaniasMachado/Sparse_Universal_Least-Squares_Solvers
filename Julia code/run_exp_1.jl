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

results_folder = "results/Experiment_$exp"

methods = ["Gurobi"]
method = methods[1]

# Gurobi parameters
opt_tol = 10^(-5)

df = DataFrame()

P1_P3_H_norm_0_list = []
P1_P3_H_norm_1_list = []
P1_P3_H_rank_list = []
P1_P3_time_list = []

PLS_H_norm_0_list = []
PLS_H_norm_1_list = []
PLS_H_rank_list = []
PLS_time_list = []

count = 0

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

    m = parse(Int, m_value)
    n = parse(Int, n_value)
    r = parse(Int, r_value)
    d = parse(Int, d_value)

    if method == "Gurobi"
        data = DataInst(A, m, n, r)

        constraints = ["P1", "P3"]

        GRB_P1_P3_time = @elapsed begin
            GRB_P1_P3_H = gurobi_solver(data, constraints, opt_tol)
        end
        GRB_P1_P3_H_norm_0 = matrix_norm_0(GRB_P1_P3_H)
        GRB_P1_P3_H_norm_1 = norm(GRB_P1_P3_H, 1)
        GRB_P1_P3_H_rank = calculate_rank(GRB_P1_P3_H)

        push!(P1_P3_H_norm_0_list, GRB_P1_P3_H_norm_0)
        push!(P1_P3_H_norm_1_list, GRB_P1_P3_H_norm_1)
        push!(P1_P3_H_rank_list, GRB_P1_P3_H_rank)
        push!(P1_P3_time_list, GRB_P1_P3_time)

        constraints = ["PLS"]

        GRB_PLS_time = @elapsed begin
            GRB_PLS_H = gurobi_solver(data, constraints, opt_tol)
        end
        GRB_PLS_H_norm_0 = matrix_norm_0(GRB_PLS_H)
        GRB_PLS_H_norm_1 = norm(GRB_PLS_H, 1)
        GRB_PLS_H_rank = calculate_rank(GRB_PLS_H)

        push!(PLS_H_norm_0_list, GRB_PLS_H_norm_0)
        push!(PLS_H_norm_1_list, GRB_PLS_H_norm_1)
        push!(PLS_H_rank_list, GRB_PLS_H_rank)
        push!(PLS_time_list, GRB_PLS_time)

        if count % 5 == 0
            GRB_P1_P3_H_norm_0_mean = mean(P1_P3_H_norm_0_list)
            GRB_P1_P3_H_norm_1_mean = mean(P1_P3_H_norm_1_list)
            GRB_P1_P3_H_rank_mean = round(Int, mean(P1_P3_H_rank_list))
            GRB_P1_P3_time_mean = mean(P1_P3_time_list)

            GRB_PLS_H_norm_0_mean = mean(PLS_H_norm_0_list)
            GRB_PLS_H_norm_1_mean = mean(PLS_H_norm_1_list)
            GRB_PLS_H_rank_mean = round(Int, mean(PLS_H_rank_list))
            GRB_PLS_time_mean = mean(PLS_time_list)

            result = DataFrame(
                m = [m],
                n = [n],
                r = [r],
                d = [d],
                GRB_P1_P3_H_norm_0_mean = [GRB_P1_P3_H_norm_0_mean],
                GRB_P1_P3_H_norm_1_mean = [GRB_P1_P3_H_norm_1_mean],
                GRB_P1_P3_H_rank_mean = [GRB_P1_P3_H_rank_mean],
                GRB_P1_P3_time_mean = [GRB_P1_P3_time_mean],
                GRB_PLS_H_norm_0_mean = [GRB_PLS_H_norm_0_mean],
                GRB_PLS_H_norm_1_mean = [GRB_PLS_H_norm_1_mean],
                GRB_PLS_H_rank_mean = [GRB_PLS_H_rank_mean],
                GRB_PLS_time_mean = [GRB_PLS_time_mean]
            )

            append!(df, result)

            empty!(P1_P3_H_norm_0_list)
            empty!(P1_P3_H_norm_1_list)
            empty!(P1_P3_H_rank_list)
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