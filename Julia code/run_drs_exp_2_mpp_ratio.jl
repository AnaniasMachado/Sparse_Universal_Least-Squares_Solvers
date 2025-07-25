using CSV
using DataFrames
using MAT
using Statistics

include("utility.jl")

exp = "2"

matrices_folder = "./Experiment_Matrices/DRS_Experiment_$(exp)"
mat_files_exp = "1"

methods = ["DRS", "DRS_Boyd", "DRS_FP"]
method = methods[3]

problems = ["P13"]
problem = problems[1]

solutions_folder = "./Solutions/DRS_Experiment_$(exp)/$(method)"

results_folder = "results/DRS_Experiment_$exp"

m_values = [60+i*20 for i in 0:27]
d_values = [100]

AMP_norm_0_div_mn_list = []

for m in m_values
    n = Int(m / 2)
    r = Int(m / 4)
    for d in d_values
        for idx in 1:5
            mat_name = "experiment_$(mat_files_exp)_matrix_m$(m)_n$(n)_r$(r)_d$(d)_idx$(idx).mat"
            mat_path = joinpath(matrices_folder, mat_name)
            mat_data = matread(mat_path)
            A = mat_data["matrix"]
            A = Matrix(A)
            AMP = pinv(A)

            AMP_norm_0_div_mn = matrix_norm_0(AMP) / (m * n)

            push!(AMP_norm_0_div_mn_list, AMP_norm_0_div_mn)
        end
    end
end

AMP_norm_0_div_mn_min = minimum(AMP_norm_0_div_mn_list)
AMP_norm_0_div_mn_mean = sum(AMP_norm_0_div_mn_list) / length(AMP_norm_0_div_mn_list)
AMP_norm_0_div_mn_max = maximum(AMP_norm_0_div_mn_list)

println("min ratio: $(AMP_norm_0_div_mn_min)")
println("mean ratio: $(AMP_norm_0_div_mn_mean)")
println("max ratio: $(AMP_norm_0_div_mn_max)")