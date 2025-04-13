using CSV
using DataFrames
using MAT
using Statistics

include("utility.jl")

exp = "6"

matrices_folder = "./Experiment_Matrices/Experiment_$(exp)"
mat_files = readdir(matrices_folder)

solvers = ["ADMM"]
solver = solvers[1]

solutions_folder = "./Solutions/Experiment_$(exp)/$(solver)"
sol_files = readdir(solutions_folder)

results_folder = "results/Experiment_$exp"

n_values = [40, 60]

solutions = Dict()

for n in n_values
    solutions[n] = []
end

# Stores the solution files of each value of m
for sol_file in sol_files
    regex = r"n_(\d+)"
    match_list = match(regex, sol_file)
    n = parse(Int, match_list.captures[1])
    push!(solutions[n], sol_file)
end

mat_files_grouped = Dict()

# Creates a list to store the experiment matrice files by value of m
for n in n_values
    mat_files_grouped[n] = []
end

# Stores the experiment matrice files by value of m
for mat_file in mat_files
    regex = r"n(\d+)"
    match_list = match(regex, mat_file)
    n = parse(Int, match_list.captures[1])
    push!(mat_files_grouped[n], mat_file)
end

stats = Dict()
stats_list = ["norm_0", "norm_1", "rank", "time"]

for n in n_values
    stats[n] = Dict()
    for stat in stats_list
        stats[n][stat] = 0
    end
end

for n in n_values
    if length(solutions[n]) < 5
        stats[n]["norm_0"] = -1.0
        stats[n]["norm_1"] = -1.0
        stats[n]["rank"] = -1.0
        stats[n]["time"] = -1.0
        continue
    end
    H_div_AMP_norm_0_sum = 0
    H_div_AMP_norm_1_sum = 0
    H_div_AMP_rank_sum = 0
    H_time_sum = 0
    for i in 1:5
        mat_file = mat_files_grouped[n][i]
        mat_path = joinpath(matrices_folder, mat_file)
        mat_data = matread(mat_path)
        A = mat_data["matrix"]
        A = Matrix(A)
        AMP = pinv(A)
        AMP_norm_0 = matrix_norm_0(AMP)
        AMP_norm_1 = norm(AMP, 1)
        AMP_rank = calculate_rank(AMP)

        sol_file = solutions[n][i]
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
    stats[n]["norm_0"] = H_div_AMP_norm_0_sum / 5
    stats[n]["norm_1"] = H_div_AMP_norm_1_sum / 5
    stats[n]["rank"] = H_div_AMP_rank_sum / 5
    stats[n]["time"] = H_time_sum / 5
end

rows = Dict{Int, Dict{String, Any}}()  # Usamos um dicionário para agrupar por instância

for (n, stats_vals) in stats
    # Verifica se a instância já existe no dicionário
    if !haskey(rows, n)
        # Se não existir, cria uma nova entrada
        rows[n] = Dict("m" => 120, "n" => n, "r" => round(n * 3/4))
    end
    
    # Adiciona os stats do problema atual na linha correspondente
    for (stat_name, stat_val) in stats_vals
        rows[n]["$(stat_name)"] = stat_val
    end
end

df = DataFrame(values(rows))

# println(names(df))

# Order dataframe by the column n
sort!(df, :n)

# Creates a list of columns from the lists problem and stats_list
dynamic_columns = [Symbol(stat) for stat in stats_list]

# Adds the fixed columns
fixed_columns = [:m, :n, :r]
all_columns = vcat(fixed_columns, dynamic_columns)

# Reorders dataframe
df = select(df, all_columns)

# Saves dataframe as a csv file
results_filename = "results_$(exp)_$(solver)_table.csv"
results_filepath = joinpath(results_folder, results_filename)
CSV.write(results_filepath, df)