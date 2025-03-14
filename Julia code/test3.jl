using CSV
using DataFrames
using MAT
using Base.GC

include("types.jl")
include("utility.jl")
include("solvers.jl")
include("drs.jl")
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
    fixed_tol = false

    DRS_time = @elapsed begin
        DRS_H = drs(A, lambda, problem, eps_opt, fixed_tol)
    end
    DRS_H_norm_0 = matrix_norm_0(DRS_H)
    DRS_H_norm_1 = norm(DRS_H, 1)

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
        DRS_time = [DRS_time]
    )

    append!(df, result)

    GC.gc()
end

CSV.write("results_drs_fixed_tol.csv", df)