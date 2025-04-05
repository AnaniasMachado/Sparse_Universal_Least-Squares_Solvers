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

epsilon = 10^(-5)

df = DataFrame()

P1_CN_list = []
PLS_CN_list = []
PMN_CN_list = []

count = 0

for mat_file in mat_files
    global count += 1

    mat_path = joinpath(matrices_folder, mat_file)
    mat_data = matread(mat_path)

    println("Computing stats of matrix: $mat_path")
    
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

    PLS_coeff = A' * A
    PMN_coeff = A * A'

    P1_S = vec(svd(A).S)
    PLS_S = vec(svd(PLS_coeff).S)
    PMN_S = vec(svd(PMN_coeff).S)

    P1_Max_SV = maximum(P1_S)
    P1_Min_SV = minimum([SV for SV in P1_S if SV > epsilon])
    PLS_Max_SV = maximum(PLS_S)
    PLS_Min_SV = minimum([SV for SV in PLS_S if SV > epsilon])
    PMN_Max_SV = maximum(PMN_S)
    PMN_Min_SV = minimum([SV for SV in PMN_S if SV > epsilon])

    P1_CN = P1_Max_SV / P1_Min_SV
    PLS_CN = PLS_Max_SV / PLS_Min_SV
    PMN_CN = PMN_Max_SV / PMN_Min_SV

    push!(P1_CN_list, P1_CN)
    push!(PLS_CN_list, PLS_CN)
    push!(PMN_CN_list, PMN_CN)

    if count % 5 == 0
        P1_CN_mean = mean(P1_CN_list)
        PLS_CN_mean = mean(PLS_CN_list)
        PMN_CN_mean = mean(PMN_CN_list)

        P1_coeff_rank = r^2
        PLS_coeff_rank = m*r
        PMN_coeff_rank = m*r

        result = DataFrame(
            m = [m],
            n = [n],
            r = [r],
            P1_coeff_rank = [P1_coeff_rank],
            PLS_coeff_rank = [PLS_coeff_rank],
            PMN_coeff_rank = [PMN_coeff_rank],
            P1_CN_mean = [P1_CN_mean],
            PLS_CN_mean = [GRB_PLS_H_norm_1_mean],
            PMN_CN_mean = [PMN_CN_mean]
        )

        append!(df, result)

        empty!(P1_CN_list)
        empty!(PLS_CN_list)
        empty!(PMN_CN_list)
    end

    GC.gc()
end

if method == "Gurobi"
    results_filename = "results_$(exp)_stats.csv"
    results_filepath = joinpath(results_folder, results_filename)
    CSV.write(results_filepath, df)
else
    throw(ErrorException("Invalid method chose."))
end