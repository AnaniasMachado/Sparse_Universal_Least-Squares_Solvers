using CSV
using DataFrames

include("types.jl")
include("utility.jl")
include("solvers.jl")
include("admm.jl")

exp = "6"

matrices_folder = "./Experiment_Matrices/Experiment_" * exp
mat_files = readdir(matrices_folder)

methods = ["Gurobi", "LS", "ADMM"]
method = [0]

# ADMM parameters
rho = 3.0
epsilon = 10^(-5)
eps_abs = epsilon
eps_rel = epsilon
fixed_tol = true
eps_opt = 10^(-5)
time_limit = 2*60*60

df = DataFrame()

for mat_file in mat_files
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
    d = parse(Int, d_value)

    A_norm_0 = matrix_norm_0(A)
    A_norm_1 = norm(A, 1)

    AMP = pinv(A)
    AMP_norm_0 = matrix_norm_0(AMP)
    AMP_norm_1 = norm(AMP, 1)

    if method == "Gurobi"
        data = DataInst(A, m, n, r)

        constraints = ["P1"]
        opt_tol = 10^(-8)

        GRB_P1_time = @elapsed begin
            GRB_P1_H = gurobi_solver(data, constraints, opt_tol)
        end
        GRB_P1_H_norm_0 = matrix_norm_0(GRB_P1_H)
        GRB_P1_H_norm_1 = norm(GRB_P1_H, 1)

        constraints = ["P1", "Sym"]
        opt_tol = 10^(-8)

        GRB_P1_Sym_time = @elapsed begin
            GRB_P1_Sym_H = gurobi_solver(data, constraints, opt_tol)
        end
        GRB_P1_Sym_H_norm_0 = matrix_norm_0(GRB_P1_Sym_H)
        GRB_P1_Sym_H_norm_1 = norm(GRB_P1_Sym_H, 1)

        result = DataFrame(
            m = [m],
            n = [n],
            r = [r],
            d = [d],
            A_norm_0 = [A_norm_0],
            A_norm_1 = [A_norm_1],
            AMP_norm_0 = [AMP_norm_0],
            AMP_norm_1 = [AMP_norm_1],
            GRB_P1_H_norm_0 = [GRB_P1_H_norm_0],
            GRB_P1_H_norm_1 = [GRB_P1_H_norm_1],
            GRB_P1_time = [GRB_P1_time],
            GRB_P1_Sym_H_norm_0 = [GRB_P1_Sym_H_norm_0],
            GRB_P1_Sym_H_norm_1 = [GRB_P1_Sym_H_norm_1],
            GRB_P1_Sym_time = [GRB_P1_Sym_time]
        )

        append!(df, result)
    elseif method == "LS"
        println("0")
    elseif method == "ADMM"
        ADMM_time = @elapsed begin
            ADMM_H = admm_p123(A, rho, eps_abs, eps_rel, fixed_tol, eps_opt, time_limit)
        end
        ADMM_H_norm_0 = matrix_norm_0(ADMM_H)
        ADMM_H_norm_1 = norm(ADMM_H, 1)

        result = DataFrame(
            m = [m],
            n = [n],
            r = [r],
            d = [d],
            A_norm_0 = [A_norm_0],
            A_norm_1 = [A_norm_1],
            AMP_norm_0 = [AMP_norm_0],
            AMP_norm_1 = [AMP_norm_1],
            ADMM_H_norm_0 = [ADMM_H_norm_0],
            ADMM_H_norm_1 = [ADMM_H_norm_1],
            ADMM_time = [ADMM_time]
        )

        append!(df, result)
    else
        throw(ErrorException("Invalid method chose."))
    end
end

if method == "Gurobi"
    CSV.write("results_GRB.csv", df)
elseif method == "LS"
    CSV.write("results_LS.csv", df)
elseif method == "ADMM"
    CSV.write("results_ADMM.csv", df)
else
    throw(ErrorException("Invalid method chose."))