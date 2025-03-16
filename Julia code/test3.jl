using CSV
using DataFrames
using MAT
using Base.GC

include("types.jl")
include("utility.jl")
include("solvers.jl")
include("drs.jl")
include("a2_basic_vec.jl")
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

    # DRS_time = @elapsed begin
    #     DRS_H = drs(A, lambda, problem, eps_opt, fixed_tol)
    # end
    DRS_time = @elapsed begin
        DRS_H, k = a2drs_basic(A, lambda, problem, eps_opt)
    end
    DRS_H_norm_0 = matrix_norm_0(DRS_H)
    DRS_H_norm_1 = norm(DRS_H, 1)

    if problem == "P123"
        println("PLS violation: $(norm(A' * A * DRS_H - A'))")
        println("P123 violation: $(norm(DRS_H * A * pinv(A) - DRS_H))")
        # if (norm(A' * A * DRS_H - A') > 10^(-5)) || (norm(DRS_H * A * pinv(A) - DRS_H) > 10^(-5))
        #     throw(ErrorException("Failed Infeasibility Check."))
        # end
    end

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

CSV.write("results_a2drs_basic.csv", df)