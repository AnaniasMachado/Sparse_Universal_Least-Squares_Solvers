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
d_values = [10, 25]

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
        solutions[problem][n] = Dict()
        for d in d_values
            solutions[problem][n][d] = []
        end
    end
end

# Stores the solution files of each value of m for each problem
for problem in problems
    for sol_file in solution_files[problem]
        regex = r"n_(\d+)"
        match_list = match(regex, sol_file)
        n = parse(Int, match_list.captures[1])

        regex = r"d_(\d+)"
        match_list = match(regex, sol_file)
        d = parse(Int, match_list.captures[1])
        push!(solutions[problem][n][d], sol_file)
    end
end

# Deletes the list of solution files that for a given value of m is less then 5
for problem in problems
    for n in n_values
        for d in d_values
            if length(solutions[problem][n][d]) < 5
                delete!(solutions[problem][n], d)
            end
        end
    end
end

mat_files_grouped = Dict()

# Creates a list to store the experiment matrice files by value of m
for n in n_values
    mat_files_grouped[n] = Dict()
    for d in d_values
        mat_files_grouped[n][d] = []
    end
end

# Stores the experiment matrice files by value of m
for mat_file in mat_files
    regex = r"n(\d+)"
    match_list = match(regex, mat_file)
    n = parse(Int, match_list.captures[1])

    regex = r"d(\d+)"
    match_list = match(regex, mat_file)
    d = parse(Int, match_list.captures[1])
    push!(mat_files_grouped[n][d], mat_file)
end

stats = Dict()
stats_list = ["A_norm_0", "AMP_norm_0", "AMP_norm_1", "hatAMP_norm_0", "hatAMP_norm_1", "norm_0", "norm_1", "rank", "time"]

for problem in problems
    stats[problem] = Dict()
    for n in n_values
        stats[problem][n] = Dict()
        for d in d_values
            stats[problem][n][d] = Dict()
            for stat in stats_list
                stats[problem][n][d][stat] = 0
            end
        end
    end
end

for problem in problems
    for n in n_values
        for d in d_values
            if !haskey(solutions[problem][n], d)
                for stat in stats_list
                    stats[problem][n][d][stat] = -1.0
                end
                continue
            end
            A_norm_0_sum = 0
            AMP_norm_0_sum = 0
            AMP_norm_1_sum = 0
            hatAMP_norm_0_sum = 0
            hatAMP_norm_1_sum = 0

            H_norm_0_sum = 0
            H_norm_1_sum = 0
            H_rank_sum = 0
            H_time_sum = 0
            for i in 1:5
                mat_file = mat_files_grouped[n][d][i]
                mat_path = joinpath(matrices_folder, mat_file)
                mat_data = matread(mat_path)
                A = mat_data["matrix"]
                A = Matrix(A)
                A_norm_0 = matrix_norm_0(A)
                AMP = pinv(A)
                AMP_norm_0 = matrix_norm_0(AMP)
                AMP_norm_1 = norm(AMP, 1)
                AMP_rank = calculate_rank(AMP)
                hatA = A' * A
                hatAMP = pinv(hatA)
                hatAMP_norm_0 = matrix_norm_0(hatAMP)
                hatAMP_norm_1 = norm(hatAMP, 1)

                sol_file = solutions[problem][n][d][i]
                sol_path = joinpath(solutions_folder, sol_file)
                sol_data = matread(sol_path)
                H = sol_data["H"]
                H = Matrix(H)
                H_norm_0 = matrix_norm_0(H)
                H_norm_1 = norm(H, 1)
                H_rank = calculate_rank(H)
                H_time = sol_data["time"]

                A_norm_0_sum += A_norm_0
                AMP_norm_0_sum += AMP_norm_0
                AMP_norm_1_sum += AMP_norm_1
                hatAMP_norm_0_sum += hatAMP_norm_0
                hatAMP_norm_1_sum += hatAMP_norm_1

                H_norm_0_sum += H_norm_0
                H_norm_1_sum += H_norm_1
                H_rank_sum += H_rank
                H_time_sum += H_time
            end
            stats[problem][n][d]["A_norm_0"] = A_norm_0_sum / 5
            stats[problem][n][d]["AMP_norm_0"] = AMP_norm_0_sum / 5
            stats[problem][n][d]["AMP_norm_1"] = AMP_norm_1_sum / 5
            stats[problem][n][d]["hatAMP_norm_0"] = hatAMP_norm_0_sum / 5
            stats[problem][n][d]["hatAMP_norm_1"] = hatAMP_norm_1_sum / 5

            stats[problem][n][d]["norm_0"] = H_norm_0_sum / 5
            stats[problem][n][d]["norm_1"] = H_norm_1_sum / 5
            stats[problem][n][d]["rank"] = H_rank_sum / 5
            stats[problem][n][d]["time"] = H_time_sum / 5
        end
    end
end

rows = Dict{Int, Dict{Int, Dict{String, Any}}}()  # Usamos um dicionário para agrupar por instância

for (problem, instances) in stats
    for (n, density_maps) in instances
        for (d, stats_vals) in density_maps
            # Verifica se a instância já existe no dicionário
            if !haskey(rows, n)
                # Se não existir, cria uma nova entrada
                rows[n] = Dict()
            end
            if !haskey(rows[n], d)
                # Se não existir, cria uma nova entrada
                rows[n][d] = Dict("m" => 120, "n" => n, "r" => round(n * 3/4), "d" => d)
            end
            
            # Adiciona os stats do problema atual na linha correspondente
            for (stat_name, stat_val) in stats_vals
                rows[n][d]["$(problem)_$(stat_name)"] = stat_val
            end
        end
    end
end

# Transformar o dicionário aninhado em uma lista de dicionários para o DataFrame
data = []
for (n, density_dict) in rows
    for (d, row_dict) in density_dict
        push!(data, row_dict)
    end
end

# Criar o DataFrame a partir da lista de dicionários
df = DataFrame(data)

# println(names(df))

# Order dataframe by the column n
sort!(df, [:n, :d])

# Creates a list of columns from the lists problem and stats_list
dynamic_columns = [Symbol(problem * "_" * stat) for problem in problems for stat in stats_list]

# Adds the fixed columns
fixed_columns = [:m, :n, :r, :d]
all_columns = vcat(fixed_columns, dynamic_columns)

# Reorders dataframe
df = select(df, all_columns)

# Saves dataframe as a csv file
results_filename = "results_$(exp)_$(solver)_table.csv"
results_filepath = joinpath(results_folder, results_filename)
CSV.write(results_filepath, df)