using CSV
using DataFrames
using MAT
using Base.GC

include("types.jl")
include("utility.jl")
include("solvers.jl")
include("drs.jl")
include("a2_basic_vec.jl")
include("a2_boyd_vec.jl")
include("a2_tr_vec.jl")
include("admm.jl")

matrix_folder = "./Experiment_Matrices/Testing_ADMM_P123_N1"
mat_files = readdir(matrix_folder)

df = DataFrame()

for mat_file in mat_files
    mat_path = joinpath(matrix_folder, mat_file)
    mat_data = matread(mat_path)
    
    A = mat_data["A"]
    A = Matrix(A)

    parts = split(mat_file, "_")

    m = parse(Int, parts[2])
    n = parse(Int, parts[3])
    r = parse(Int, split(parts[4], ".")[1])

    A_norm_0 = matrix_norm_0(A)
    A_norm_1 = norm(A, 1)

    AMP = pinv(A)
    AMP_norm_0 = matrix_norm_0(AMP)
    AMP_norm_1 = norm(AMP, 1)

    lambda = 10^(-2)
    problem = "P123"
    eps_opt = 10^(-5)
    fixed_tol = true

    DRS_time = @elapsed begin
        DRS_H, DRS_k = drs(A, lambda, problem, eps_opt)
    end
    DRS_H_norm_0 = matrix_norm_0(DRS_H)
    DRS_H_norm_1 = norm(DRS_H, 1)

    GC.gc()

    A2DRS_Basic_time = @elapsed begin
        A2DRS_Basic_H, A2DRS_Basic_k = a2drs_basic(A, lambda, problem, eps_opt)
    end
    A2DRS_Basic_H_norm_0 = matrix_norm_0(A2DRS_Basic_H)
    A2DRS_Basic_H_norm_1 = norm(A2DRS_Basic_H, 1)

    GC.gc()

    A2DRS_Boyd_time = @elapsed begin
        A2DRS_Boyd_H, A2DRS_Boyd_k = a2drs_boyd(A, lambda, problem, eps_opt)
    end
    A2DRS_Boyd_H_norm_0 = matrix_norm_0(A2DRS_Boyd_H)
    A2DRS_Boyd_H_norm_1 = norm(A2DRS_Boyd_H, 1)

    GC.gc()

    A2DRS_TR_time = @elapsed begin
        A2DRS_TR_H, A2DRS_TR_k = a2drs_tr(A, lambda, problem, eps_opt)
    end
    A2DRS_TR_H_norm_0 = matrix_norm_0(A2DRS_TR_H)
    A2DRS_TR_H_norm_1 = norm(A2DRS_TR_H, 1)

    GC.gc()

    result = DataFrame(
        m = [m],
        n = [n],
        r = [r],
        A_norm_0 = [A_norm_0],
        A_norm_1 = [A_norm_1],
        AMP_norm_0 = [AMP_norm_0],
        AMP_norm_1 = [AMP_norm_1],
        DRS_H_norm_0 = [DRS_H_norm_0],
        DRS_H_norm_1 = [DRS_H_norm_1],
        DRS_time = [DRS_time],
        DRS_k = [DRS_k],
        A2DRS_Basic_H_norm_0 = [A2DRS_Basic_H_norm_0],
        A2DRS_Basic_H_norm_1 = [A2DRS_Basic_H_norm_1],
        A2DRS_Basic_time = [A2DRS_Basic_time],
        A2DRS_Basic_k = [A2DRS_Basic_k],
        A2DRS_Boyd_H_norm_0 = [A2DRS_Boyd_H_norm_0],
        A2DRS_Boyd_H_norm_1 = [A2DRS_Boyd_H_norm_1],
        A2DRS_Boyd_time = [A2DRS_Boyd_time],
        A2DRS_Boyd_k = [A2DRS_Boyd_k],
        A2DRS_TR_H_norm_0 = [A2DRS_TR_H_norm_0],
        A2DRS_TR_H_norm_1 = [A2DRS_TR_H_norm_1],
        A2DRS_TR_time = [A2DRS_TR_time],
        A2DRS_TR_k = [A2DRS_TR_k],
    )

    append!(df, result)

    GC.gc()
end

CSV.write("results_drs_a2drs_comparison.csv", df)