using CSV
using DataFrames
using MAT
using Statistics

include("utility.jl")

exp = "1"

matrices_folder = "./Experiment_Matrices/Experiment_$(exp)"
mat_files = readdir(matrices_folder)

solutions_folder = "./Solutions/Experiment_$(exp)/Gurobi"
sol_files = readdir(solutions_folder)

results_folder = "results/Experiment_$exp"

problems = ["P1_P3", "PLS"]

m_values = [60+i*20 for i in 0:27]

# Creates a list to store the solution files of each problem
solution_files = Dict()
for problem in problems
    solution_files[problem] = []
end

# Stores the solution files of each problem in its own list
for sol_file in sol_files
    for problem in problems
        if occursin(problem, sol_file)
            push!(solution_files[problem], sol_file)
        end
    end
end

solutions = Dict()

# Creates a list to store the solution files of each value of m for each problem
for problem in problems
    solutions[problem] = Dict()
    for m in m_values
        solutions[problem][m] = []
    end
end

# Stores the solution files of each value of m for each problem
for problem in problems
    for sol_file in solution_files[problem]
        regex = r"m_(\d+)"
        match_list = match(regex, sol_file)
        m = parse(Int, match_list.captures[1])
        push!(solutions[problem][m], sol_file)
    end
end

# Deletes the list of solution files that for a given value of m is less then 5
for problem in problems
    for m in m_values
        if length(solutions[problem][m]) < 5
            delete!(solutions[problem], m)
        end
    end
end

mat_files_grouped = Dict()

# Creates a list to store the experiment matrice files by value of m
for m in m_values
    mat_files_grouped[m] = []
end

# Stores the experiment matrice files by value of m
for mat_file in mat_files
    regex = r"m(\d+)"
    match_list = match(regex, mat_file)
    m = parse(Int, match_list.captures[1])
    push!(mat_files_grouped[m], mat_file)
end

stats = Dict()
stats_list = ["norm_0", "norm_1", "rank", "time"]

for problem in problems
    stats[problem] = Dict()
    for m in m_values
        stats[problem][m] = Dict()
        for stat in stats_list
            stats[problem][m][stat] = 0
        end
    end
end

for problem in problems
    for m in m_values
        if !haskey(solutions[problem], m)
            stats[problem][m]["norm_0"] = -1.0
            stats[problem][m]["norm_1"] = -1.0
            stats[problem][m]["rank"] = -1.0
            stats[problem][m]["time"] = -1.0
            continue
        end
        H_div_AMP_norm_0_sum = 0
        H_div_AMP_norm_1_sum = 0
        H_div_AMP_rank_sum = 0
        H_time_sum = 0
        for i in 1:5
            mat_file = mat_files_grouped[m][i]
            mat_path = joinpath(matrices_folder, mat_file)
            mat_data = matread(mat_path)
            A = mat_data["matrix"]
            A = Matrix(A)
            AMP = pinv(A)
            AMP_norm_0 = matrix_norm_0(AMP)
            AMP_norm_1 = norm(AMP, 1)
            AMP_rank = calculate_rank(AMP)

            sol_file = solutions[problem][m][i]
            sol_path = joinpath(solutions_folder, sol_file)
            sol_data = matread(sol_path)
            H = sol_data["H"]
            H = Matrix(H)
            H_norm_0 = matrix_norm_0(H)
            H_norm_1 = norm(H, 1)
            H_rank = calculate_rank(H)
            H_time = sol_data["time"]

            H_div_AMP_norm_0_sum += H_norm_0 / AMP_norm_0
            H_div_AMP_norm_1_sum += H_norm_1 / AMP_norm_1
            H_div_AMP_rank_sum += H_rank / AMP_rank
            H_time_sum += H_time
        end
        stats[problem][m]["norm_0"] = H_div_AMP_norm_0_sum / 5
        stats[problem][m]["norm_1"] = H_div_AMP_norm_1_sum / 5
        stats[problem][m]["rank"] = H_div_AMP_rank_sum / 5
        stats[problem][m]["time"] = H_time_sum / 5
    end
end

rows = Dict{Int, Dict{String, Any}}()  # Usamos um dicionário para agrupar por instância

for (problem, instances) in stats
    for (m, stats_vals) in instances
        # Verifica se a instância já existe no dicionário
        if !haskey(rows, m)
            # Se não existir, cria uma nova entrada
            rows[m] = Dict("m" => m, "n" => m / 2, "r" => m / 4)
        end
        
        # Adiciona os stats do problema atual na linha correspondente
        for (stat_name, stat_val) in stats_vals
            rows[m]["$(problem)_$(stat_name)"] = stat_val
        end
    end
end

df = DataFrame(values(rows))

# println(names(df))

# Order dataframe by the column m
sort!(df, :m)

# Creates a list of columns from the lists problem and stats_list
dynamic_columns = [Symbol(problem * "_" * stat) for problem in problems for stat in stats_list]

# Adds the fixed columns
fixed_columns = [:m, :n, :r]
all_columns = vcat(fixed_columns, dynamic_columns)

# Reorders dataframe
df = select(df, all_columns)

# Saves dataframe as a csv file
results_filename = "results_$(exp)_table.csv"
results_filepath = joinpath(results_folder, results_filename)
CSV.write(results_filepath, df)