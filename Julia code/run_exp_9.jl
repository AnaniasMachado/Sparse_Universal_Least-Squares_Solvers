using CSV
using DataFrames
using MAT
using Base.GC
using Statistics

include("types.jl")
include("utility.jl")
include("solvers.jl")
include("solvers_cal.jl")

exp = "9"

matrices_folder = "./Experiment_Matrices/Experiment_" * exp
mat_files = readdir(matrices_folder)

results_folder = "results/Experiment_$exp"

solutions_folder = "./Solutions/Experiment_" * exp

methods = ["Gurobi", "Gurobi_Cal"]
method = methods[1]

# Gurobi parameter
constraints_set = [["P1", "P3"], ["P13R"], ["PLS"]]
constraints = constraints_set[2]
problems = ["P13"]
problem = problems[1]

eps_opt = 10^(-5)
time_limit = 1200

df = DataFrame()

H_div_mr_norm_0_list = []
H_div_AMP_norm_0_list = []
H_div_AMP_norm_1_list = []
H_rank_ratio_list = []
time_list = []

count = 0
num_instances = 5
min_unsolvable_m = Dict()

for mat_file in mat_files
    global count += 1

    mat_path = joinpath(matrices_folder, mat_file)
    mat_data = matread(mat_path)

    println("Solving for matrix: $mat_path")
    
    A = mat_data["matrix"]
    A = Matrix(A)
    AMP = pinv(A)

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

    data = DataInst(A, m, n, r, AMP=AMP)

    if !haskey(min_unsolvable_m, d)
        min_unsolvable_m[d] = Inf
    end

    H_div_mr_norm_0 = -1.0
    H_div_AMP_norm_0 = -1.0
    H_div_AMP_norm_1 = -1.0
    H_rank_ratio = -1.0
    time = -1.0
    if (m < min_unsolvable_m[d])
        if method == "Gurobi"
            try
                time = @elapsed begin
                    H = gurobi_solver(data, constraints, eps_opt, time_limit)
                end
                H_norm_0 = matrix_norm_0(H)
                H_norm_1 = norm(H, 1)
                H_rank = calculate_rank(H)

                H_div_mr_norm_0 = H_norm_0 / (m * r)
                H_div_AMP_norm_0 = H_norm_0 / matrix_norm_0(AMP)
                H_div_AMP_norm_1 = H_norm_1 / norm(AMP, 1)
                H_rank_ratio = H_rank / r

                problem_label = join(constraints, "_")

                solution_filename = "Gurobi/Experiment_$(exp)_$(problem_label)_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                solution_filepath = joinpath(solutions_folder, solution_filename)
                matwrite(solution_filepath, Dict("H" => H, "time" => time))
            catch e
                if isa(e, ErrorException)
                    global min_unsolvable_m[d] = min(m, min_unsolvable_m[d])
                else
                    throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
                end
            end
        elseif method == "Gurobi_Cal"
            try
                time = @elapsed begin
                    H = gurobi_solver_cal(data, problem, eps_opt, time_limit)
                end
                H_norm_0 = matrix_norm_0(H)
                H_norm_1 = norm(H, 1)
                H_rank = calculate_rank(H)

                H_div_mr_norm_0 = H_norm_0 / (m * r)
                H_div_AMP_norm_0 = H_norm_0 / matrix_norm_0(AMP)
                H_div_AMP_norm_1 = H_norm_1 / norm(AMP, 1)
                H_rank_ratio = H_rank / r

                solution_filename = "Gurobi_Cal/Experiment_$(exp)_$(problem)_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                solution_filepath = joinpath(solutions_folder, solution_filename)
                matwrite(solution_filepath, Dict("H" => H, "time" => time))
            catch e
                if isa(e, ErrorException)
                    global min_unsolvable_m[d] = min(m, min_unsolvable_m[d])
                else
                    throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
                end
            end
        else
            throw(ErrorException("Invalid method chose."))
        end
    end

    push!(H_div_mr_norm_0_list, H_div_mr_norm_0)
    push!(H_div_AMP_norm_0_list, H_div_AMP_norm_0)
    push!(H_div_AMP_norm_1_list, H_div_AMP_norm_1)
    push!(H_rank_ratio_list, H_rank_ratio)
    push!(time_list, time)

    GC.gc()

    if count % num_instances == 0
        H_div_mr_norm_0_mean = -1.0
        H_div_AMP_norm_0_mean = -1.0
        H_div_AMP_norm_1_mean = -1.0
        H_rank_ratio_mean = -1.0
        time_mean = -1.0

        if !(-1.0 in H_div_AMP_norm_0_list)
            H_div_mr_norm_0_mean = mean(H_div_mr_norm_0_list)
            H_div_AMP_norm_0_mean = mean(H_div_AMP_norm_0_list)
            H_div_AMP_norm_1_mean = mean(H_div_AMP_norm_1_list)
            H_rank_ratio_mean = mean(H_rank_ratio_list)
            time_mean = mean(time_list)
        end

        result = DataFrame(
            m = [m],
            r = [r],
            d = [d],
            A_norm_0 = [matrix_norm_0(A)],
            A_norm_1 = [norm(A, 1)],
            AMP_norm_0 = [matrix_norm_0(AMP)],
            AMP_norm_1 = [norm(AMP, 1)],
            H_div_mr_norm_0_mean = [H_div_mr_norm_0_mean],
            H_div_AMP_norm_0_mean = [H_div_AMP_norm_0_mean],
            H_div_AMP_norm_1_mean = [H_div_AMP_norm_1_mean],
            H_rank_ratio_mean = [H_rank_ratio_mean],
            time_mean = [time_mean]
        )

        append!(df, result)

        empty!(H_div_mr_norm_0_list)
        empty!(H_div_AMP_norm_0_list)
        empty!(H_div_AMP_norm_1_list)
        empty!(H_rank_ratio_list)
        empty!(time_list)

        GC.gc()
    end
end

if method == "Gurobi"
    problem_label = join(constraints, "_")
    results_filename = "results_$(exp)_Gurobi_$(problem_label).csv"
    results_filepath = joinpath(results_folder, results_filename)
    CSV.write(results_filepath, df)
elseif method == "Gurobi_Cal"
    results_filename = "results_$(exp)_Gurobi_Cal_$(problem).csv"
    results_filepath = joinpath(results_folder, results_filename)
    CSV.write(results_filepath, df)
else
    throw(ErrorException("Invalid method chose."))
end