using LinearAlgebra
using JuMP, MosekTools, Gurobi
using MAT,DelimitedFiles
using SparseArrays
# import MutableArithmetics
# const MA = MutableArithmetics
# const SA = MA.SparseArrays
# using Base.Threads

using JLD
using Plots, LaTeXStrings
using Printf

include("types.jl")
include("util.jl")
include("closed_form_sols.jl")
include("admm.jl")
include("solvers.jl")
include("local_search.jl")

arr_resADMM_1,arr_resADMM_0 = init_ginv_res_admm(),init_ginv_res_admm()
println("\nStarting procedure...")


fixed_tol = false

fixed_tol ? TP = :rel_tol : TP=:fixed_tol

expn = "6"

matrix_folder = "./Experiment_Matrices/Experiment_" * expn
mat_files = readdir(matrix_folder)

result_matrix_folder = "./ADMM_Matrices/Experiment_" * expn

for mat_file in mat_files[1:2]
    if endswith(mat_file, ".mat")
        mat_path = joinpath(matrix_folder, mat_file)
        mat_data = matread(mat_path)
    
        A = mat_data["matrix"]
        A = Matrix(A)

        m_value = match(r"m(\d+)", mat_file).captures[1]
        n_value = match(r"n(\d+)", mat_file).captures[1]
        r_value = match(r"r(\d+)", mat_file).captures[1]
        d_value = match(r"d(\d+)", mat_file).captures[1]

        m = parse(Int, m_value)
        n = parse(Int, n_value)
        r = parse(Int, r_value)
        dens_int = parse(Int, d_value)

        inst = GinvInst(A,m,n,r);
        ginvInit = getInitialInfoGinv(inst)
        if fixed_tol
            # Run ADMM based on the 1-norm
            time_admm_1 = @elapsed admmsol_1 = admm1norm(ginvInit;eps_opt=1e-2,stop_limit=:OptGap);
            result_mat_file = replace(mat_file, "matrix" => "matrixADMM")
        else
            # Run ADMM based on the 1-norm
            time_admm_1 = @elapsed admmsol_1 = admm1norm(ginvInit;eps_abs=1e-4,eps_rel=1e-4,stop_limit=:Boyd);
            result_mat_file = replace(mat_file, "matrix" => "matrixADMMe")
        end
        result_mat_path = joinpath(result_matrix_folder, result_mat_file)
        matwrite(result_mat_path, Dict("matrix" => admmsol_1.H))
    end
end