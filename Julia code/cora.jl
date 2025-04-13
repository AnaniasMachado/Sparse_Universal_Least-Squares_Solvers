include("types.jl")
include("utility.jl")
include("solvers.jl")

using CSV
using DataFrames
using Base.GC

file_path = "./Experiment_Matrices/cora/cora.content"

data = CSV.File(file_path, delim='\t', header=false)
df = DataFrame(data)

filtered_df = df[:, 2:end-1]

A = Matrix(filtered_df)

m, n = size(A)
r = calculate_rank(A)
d = matrix_norm_0(A) / (m * n)

println("m = $m")
println("n = $n")
println("r = $r")
println("d = $d")

function vector_norm_0(x)
    norm_0 = 0
    for i in 1:length(x)
        if x[i] != 0
            norm_0 += 1
        end
    end
    return norm_0
end

function k_sparsest_columns(A, k)
    sparsity = [vector_norm_0(A[:, j]) for j in 1:n]
    
    sorted_indices = sortperm(sparsity)[1:k]
    
    return A[:, sorted_indices[1:k]]
end

function remove_k_densest_columns(A, k)
    sparsity = [vector_norm_0(A[:, j]) for j in 1:n]
    
    sorted_indices = sortperm(sparsity, rev=true)[1:k]
    
    return A[:, Not(sorted_indices[1:k])]
end

k = 800
E = remove_k_densest_columns(A, k)

m, n = size(E)
r = calculate_rank(E)
d = matrix_norm_0(E) / (m * n)

println("m = $m")
println("n = $n")
println("r = $r")
println("d = $d")

k = 600
B = k_sparsest_columns(E, k)
p = round(Int, k / 2)
C = zeros(m, p)
for i in 1:p
    C[:, i] = B[:, i*2-1] + B[:, i*2]
end

D = hcat(E, C)

m, n = size(D)
r = calculate_rank(D)
d = matrix_norm_0(D) / (m * n)

k = 600
F = k_sparsest_columns(D, k)
p = round(Int, k / 2)
G = zeros(m, p)
for i in 1:p
    G[:, i] = F[:, i*2-1] + F[:, i*2]
end

K = hcat(D, G)

m, n = size(K)
r = calculate_rank(K)
d = matrix_norm_0(K) / (m * n)

println("m = $m")
println("n = $n")
println("r = $r")
println("d = $d")

results_folder = "results/cora"

methods = ["Gurobi"]
method = methods[1]

# Gurobi parameters
opt_tol = 10^(-5)
time_limit = 12*60*60

df = DataFrame()

if method == "Gurobi"
    data = DataInst(K, m, n, r)

    constraints = ["P1", "P3"]

    GRB_P1_P3_time = -1.0
    try
        global GRB_P1_P3_time = @elapsed begin
            GRB_P1_P3_H = gurobi_solver(data, constraints, opt_tol, time_limit)
        end
    catch e
        if isa(e, ErrorException)
            global GRB_P1_P3_time = -1.0
        else
            throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
        end
    end

    GC.gc()

    constraints = ["PLS"]

    GRB_PLS_time = -1.0
    GRB_PLS_H_norm_0 = -1.0
    GRB_PLS_H_norm_1 = -1.0
    GRB_PLS_H_rank = -1.0

    try
        global GRB_PLS_time = @elapsed begin
            GRB_PLS_H = gurobi_solver(data, constraints, opt_tol, time_limit)
        end
        global GRB_PLS_H_norm_0 = matrix_norm_0(GRB_PLS_H)
        global GRB_PLS_H_norm_1 = norm(GRB_PLS_H, 1)
        global GRB_PLS_H_rank = calculate_rank(GRB_PLS_H)
    catch e
        if isa(e, ErrorException)
            global GRB_PLS_time = -1.0
        else
            throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
        end
    end

    result = DataFrame(
        m = [m],
        n = [n],
        r = [r],
        d = [d],
        GRB_PLS_H_norm_0 = [GRB_PLS_H_norm_0],
        GRB_PLS_H_norm_1 = [GRB_PLS_H_norm_1],
        GRB_PLS_H_rank = [GRB_PLS_H_rank],
        GRB_P1_P3_time_mean = [GRB_P1_P3_time_mean],
        GRB_PLS_time = [GRB_PLS_time]
    )

    append!(df, result)

    GC.gc()
else
    throw(ErrorException("Invalid method chose."))
    end
end

if method == "Gurobi"
    results_filename = "results_cora_GRB.csv"
    results_filepath = joinpath(results_folder, results_filename)
    CSV.write(results_filepath, df)
else
    throw(ErrorException("Invalid method chose."))
end