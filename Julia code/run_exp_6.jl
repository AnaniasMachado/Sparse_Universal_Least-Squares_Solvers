using CSV
using DataFrames
using MAT
using Base.GC

include("types.jl")
include("utility.jl")
include("solvers.jl")
include("admm.jl")

exp = "6"

matrices_folder = "./Experiment_Matrices/Experiment_" * exp
mat_files = readdir(matrices_folder)

solutions_folder = "./Solutions/Experiment_" * exp

results_folder = "results/Experiment_$exp"

methods = ["Gurobi", "LS", "ADMM", "hatAMP_data"]
method = methods[3]

# Gurobi parameters
opt_tol = 10^(-5)

# ADMM parameters
rho = 3.0
epsilon = 10^(-5)
eps_abs = epsilon
eps_rel = epsilon
fixed_tol = false
eps_opt = epsilon
time_limit = 1200

df = DataFrame()

# for mat_file in mat_files[41:end]
for mat_file in mat_files
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

    A_norm_0 = matrix_norm_0(A)
    A_norm_1 = norm(A, 1)

    AMP = pinv(A)
    AMP_norm_0 = matrix_norm_0(AMP)
    AMP_norm_1 = norm(AMP, 1)

    if method == "Gurobi"
        hatA = A' * A
        hatm, hatn = size(hatA)
        data = DataInst(hatA, hatm, hatn, r)

        # constraints = ["P1"]

        # GRB_P1_time = @elapsed begin
        #     GRB_P1_H = gurobi_solver(data, constraints, opt_tol)
        # end
        # GRB_P1_H_norm_0 = matrix_norm_0(GRB_P1_H)
        # GRB_P1_H_norm_1 = norm(GRB_P1_H, 1)
        # solution_filename = "Gurobi/Experiment_$(exp)_P1_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
        # solution_filepath = joinpath(solutions_folder, solution_filename)
        # matwrite(solution_filepath, Dict("H" => GRB_P1_H, "time" => GRB_P1_time))

        constraints = ["P1", "Sym"]

        GRB_P1_Sym_time = @elapsed begin
            GRB_P1_Sym_H = gurobi_solver(data, constraints, opt_tol)
        end
        GRB_P1_Sym_H_norm_0 = matrix_norm_0(GRB_P1_Sym_H)
        GRB_P1_Sym_H_norm_1 = norm(GRB_P1_Sym_H, 1)
        solution_filename = "Gurobi/Experiment_$(exp)_P1_Sym_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
        solution_filepath = joinpath(solutions_folder, solution_filename)
        matwrite(solution_filepath, Dict("H" => GRB_P1_Sym_H, "time" => GRB_P1_Sym_time))

        # result = DataFrame(
        #     m = [m],
        #     n = [n],
        #     r = [r],
        #     d = [d],
        #     A_norm_0 = [A_norm_0],
        #     A_norm_1 = [A_norm_1],
        #     AMP_norm_0 = [AMP_norm_0],
        #     AMP_norm_1 = [AMP_norm_1],
        #     GRB_P1_H_norm_0 = [GRB_P1_H_norm_0],
        #     GRB_P1_H_norm_1 = [GRB_P1_H_norm_1],
        #     GRB_P1_time = [GRB_P1_time],
        #     GRB_P1_Sym_H_norm_0 = [GRB_P1_Sym_H_norm_0],
        #     GRB_P1_Sym_H_norm_1 = [GRB_P1_Sym_H_norm_1],
        #     GRB_P1_Sym_time = [GRB_P1_Sym_time]
        # )

        result = DataFrame(
            m = [m],
            n = [n],
            r = [r],
            d = [d],
            A_norm_0 = [A_norm_0],
            A_norm_1 = [A_norm_1],
            AMP_norm_0 = [AMP_norm_0],
            AMP_norm_1 = [AMP_norm_1],
            GRB_P1_Sym_H_norm_0 = [GRB_P1_Sym_H_norm_0],
            GRB_P1_Sym_H_norm_1 = [GRB_P1_Sym_H_norm_1],
            GRB_P1_Sym_time = [GRB_P1_Sym_time]
        )

        append!(df, result)

        GC.gc()
    elseif method == "LS"
        println("0")
    elseif method == "ADMM"
        ADMM_time = @elapsed begin
            ADMM_H = admm_p123(A, rho, eps_abs, eps_rel, fixed_tol, eps_opt, time_limit)
        end
        ADMM_H_norm_0 = matrix_norm_0(ADMM_H)
        ADMM_H_norm_1 = norm(ADMM_H, 1)
        if fixed_tol
            solution_filename = "ADMMe/Experiment_$(exp)_P123_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
            solution_filepath = joinpath(solutions_folder, solution_filename)
            matwrite(solution_filepath, Dict("H" => ADMM_H, "time" => ADMM_time))
        else
            solution_filename = "ADMM/Experiment_$(exp)_P123_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
            solution_filepath = joinpath(solutions_folder, solution_filename)
            matwrite(solution_filepath, Dict("H" => ADMM_H, "time" => ADMM_time))
        end

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

        GC.gc()
    elseif method == "hatAMP_data"
        hatA = A' * A
        hatAMP = pinv(hatA)

        hatAMP_norm_0 = matrix_norm_0(hatAMP)
        hatAMP_norm1 = norm(hatAMP, 1)

        result = DataFrame(
            hatAMP_norm_0 = hatAMP_norm_0,
            hatAMP_norm1 = hatAMP_norm1
        )

        append!(df, result)
    else
        throw(ErrorException("Invalid method chose."))
    end
end

if method == "Gurobi"
    results_filename = "results_$(exp)_GRB.csv"
    results_filepath = joinpath(results_folder, results_filename)
    CSV.write(results_filepath, df)
elseif method == "LS"
    results_filename = "results_$(exp)_LS.csv"
    results_filepath = joinpath(results_folder, results_filename)
    CSV.write(results_filepath, df)
elseif method == "ADMM"
    if fixed_tol
        results_filename = "results_$(exp)_ADMMe.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    else
        results_filename = "results_$(exp)_ADMM.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    end
elseif method == "hatAMP_data"
        results_filename = "results_$(exp)_hatAMP.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
else
    throw(ErrorException("Invalid method chose."))
end