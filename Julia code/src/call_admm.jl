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

arr_resADMM_1, arr_resADMM_0 = init_ginv_res_admm(), init_ginv_res_admm()

function run_admm_p123(A::Matrix{Float64}, m::Int64, n::Int64, r::Int64, rho::Float64, eps_abs::Float64, eps_rel::Float64, eps_opt::Float64, fixed_tol::Bool)
    inst = GinvInst(A,m,n,r);
    ginvInit = getInitialInfoGinv(inst)
    if fixed_tol
        # Run ADMM based on the 1-norm
        admmsol_1 = admm1norm(ginvInit;eps_opt=eps_opt,stop_limit=:OptGap);
        return admmsol_1
    else
        # Run ADMM based on the 1-norm
        admmsol_1 = admm1norm(ginvInit;eps_abs=eps_abs,eps_rel=eps_rel,stop_limit=:Boyd);
        return admmsol_1
    end
end