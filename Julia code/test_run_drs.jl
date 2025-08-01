using MAT

include("utility.jl")
include("types.jl")
include("drs.jl")

m = 100
n = 50
r = 25

A = gen_random_rank_r_matrix(m, n, r)

matrix_folder = "./Experiment_Matrices/DRS_Experiment_2"
mat_files = readdir(matrix_folder)
mat_file = mat_files[1]
mat_path = joinpath(matrix_folder, mat_file)
mat_data = matread(mat_path)
A = mat_data["matrix"]
A = Matrix(A)
AMP = pinv(A)
m, n = size(A)
r = rank(A)

println("m = $m, n = $n, r = $r")

data = DataInst(A, m, n, r, AMP=AMP)
constraints = ["PLS"]
problem = "PLS"
rho = 3.0
lambda = 10^(-2)

epsilon = 10^(-5)
eps_opt = epsilon
eps_abs = epsilon
eps_rel = 10^(-4)
fixed_tol = false

time_limit = 1200
stop_crits = ["Epsilon", "Boyd", "Fixed_Point", "FP_Epsilon", "FP_Soft"]
stop_crit = stop_crits[5]

DRS_time = @elapsed begin
    DRS_H, DRS_k = drs(A, lambda, eps_abs, eps_rel, problem, fixed_tol, eps_opt, stop_crit, time_limit)
end
DRS_H_norm_0 = matrix_norm_0(DRS_H)
DRS_H_norm_1 = norm(DRS_H, 1)

println("Stop crit: $stop_crit")
println("DRS time: $DRS_time")
println("DRS k: $DRS_k")
println("DRS norm 1: $DRS_H_norm_1")
println("DRS norm 0: $DRS_H_norm_0")
println("DRS bound ratio: $(DRS_H_norm_0 / (m * r))")