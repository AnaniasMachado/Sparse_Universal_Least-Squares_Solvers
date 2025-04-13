using CSV
using DataFrames
using MAT
using Statistics

include("utility.jl")

exp = "6"

matrices_folder = "./Experiment_Matrices/Experiment_$(exp)"
mat_files = readdir(matrices_folder)

solvers = ["LS"]
solver = solvers[1]

solutions_folder = "./Solutions/Experiment_$(exp)/$(solver)"
sol_files = readdir(solutions_folder)

results_folder = "results/Experiment_$exp"

problems = ["P123", "P1_Sym"]

n_values = [40, 60]

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
    for n in n_values
        solutions[problem][n] = []
    end
end

# Stores the solution files of each value of m for each problem
for problem in problems
    for sol_file in solution_files[problem]
        regex = r"n_(\d+)"
        match_list = match(regex, sol_file)
        n = parse(Int, match_list.captures[1])
        push!(solutions[problem][n], sol_file)
    end
end

# Deletes the list of solution files that for a given value of m is less then 5
for problem in problems
    for n in n_values
        if length(solutions[problem][n]) < 5
            delete!(solutions[problem], n)
        end
    end
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

for problem in problems
    stats[problem] = Dict()
    for n in n_values
        stats[problem][n] = Dict()
        for stat in stats_list
            stats[problem][n][stat] = 0
        end
    end
end

for problem in problems
    for n in n_values
        if !haskey(solutions[problem], n)
            stats[problem][n]["norm_0"] = -1.0
            stats[problem][n]["norm_1"] = -1.0
            stats[problem][n]["rank"] = -1.0
            stats[problem][n]["time"] = -1.0
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

            sol_file = solutions[problem][n][i]
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
        stats[problem][n]["norm_0"] = H_div_AMP_norm_0_sum / 5
        stats[problem][n]["norm_1"] = H_div_AMP_norm_1_sum / 5
        stats[problem][n]["rank"] = H_div_AMP_rank_sum / 5
        stats[problem][n]["time"] = H_time_sum / 5
    end
end

rows = Dict{Int, Dict{String, Any}}()  # Usamos um dicionário para agrupar por instância

for (problem, instances) in stats
    for (n, stats_vals) in instances
        # Verifica se a instância já existe no dicionário
        if !haskey(rows, n)
            # Se não existir, cria uma nova entrada
            rows[n] = Dict("m" => 120, "n" => n, "r" => round(n * 3/4))
        end
        
        # Adiciona os stats do problema atual na linha correspondente
        for (stat_name, stat_val) in stats_vals
            rows[n]["$(problem)_$(stat_name)"] = stat_val
        end
    end
end

df = DataFrame(values(rows))

# println(names(df))

# Order dataframe by the column n
sort!(df, :n)

# Creates a list of columns from the lists problem and stats_list
dynamic_columns = [Symbol(problem * "_" * stat) for problem in problems for stat in stats_list]

# Adds the fixed columns
fixed_columns = [:m, :n, :r]
all_columns = vcat(fixed_columns, dynamic_columns)

# Reorders dataframe
df = select(df, all_columns)

# Saves dataframe as a csv file
results_filename = "results_$(exp)_$(solver)_table.csv"
results_filepath = joinpath(results_folder, results_filename)
CSV.write(results_filepath, df)