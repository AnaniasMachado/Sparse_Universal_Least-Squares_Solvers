using MAT

include("utility.jl")
include("types.jl")
include("solvers.jl")
include("drs.jl")
include("a2_boyd.jl")
include("a2_basic.jl")
include("admm.jl")

# m = 10
# n = 5
# r = 2

# m = 50
# n = 25
# r = 10

m = 100
n = 50
r = 25

# A = gen_random_rank_r_matrix(m, n, r)

matrix_folder = "./Experiment_Matrices/Testing_ADMM_P123_N1"
mat_files = readdir(matrix_folder)
mat_file = mat_files[1]
mat_path = joinpath(matrix_folder, mat_file)
mat_data = matread(mat_path)
A = mat_data["A"]
A = Matrix(A)

# mat_path = "matrixA.mat"
# mat_data = matread(mat_path)
# A = mat_data["A"]
# A = Matrix(A)

data = DataInst(A, m, n, r)
constraints = ["PLS"]
problem = "P123"
eps_opt = 10^(-5)
# lambda = 0.28
# lambda = 0.000026
# lambda = 0.0026
lambda = 10^(-2)

# lambda = e-4; time = 38,9s; quality = slightly worse than optimal
# lambda = e-3; time = 51,1s; quality = slightly better than optimal
# lambda = e-2; time = 47,2s; quality = slightly better than optimal

GRB_time = @elapsed begin
    GRB_H = gurobi_solver(data, constraints, eps_opt)
end
GRB_H_norm_0 = matrix_norm_0(GRB_H)
GRB_H_norm_1 = norm(GRB_H, 1)

DRS_time = @elapsed begin
    DRS_H = drs(A, lambda, problem, eps_opt)
end
DRS_H_norm_0 = matrix_norm_0(DRS_H)
DRS_H_norm_1 = norm(DRS_H, 1)

if problem == "PLS"
    if norm(A' * A * DRS_H - A') > 10^(-5)
        throw(ErrorException("Failed Infeasibility Check."))
    end
elseif problem == "P123"
    println("PLS violation: $(norm(A' * A * DRS_H - A'))")
    println("P123 violation: $(norm(DRS_H * A * pinv(A) - DRS_H))")
    if norm(A' * A * DRS_H - A') + norm(DRS_H * A * pinv(A) - DRS_H) > 10^(-5)
        throw(ErrorException("Failed Infeasibility Check."))
    end
end

# A2DRS_time = @elapsed begin
#     A2DRS_H = a2drs_tr(A, lambda, problem, eps_opt)
# end
# A2DRS_H_norm_0 = matrix_norm_0(A2DRS_H)
# A2DRS_H_norm_1 = norm(A2DRS_H, 1)

A2DRS_Boyd_time = @elapsed begin
    A2DRS_Boyd_H = a2drs_boyd(A, lambda, problem, eps_opt)
end
A2DRS_Boyd_H_norm_0 = matrix_norm_0(A2DRS_Boyd_H)
A2DRS_Boyd_H_norm_1 = norm(A2DRS_Boyd_H, 1)

# A2DRS_ATM_time = @elapsed begin
#     A2DRS_ATM_H = a2drs_atm(A, lambda, problem, eps_opt)
# end
# A2DRS_ATM_H_norm_0 = matrix_norm_0(A2DRS_ATM_H)
# A2DRS_ATM_H_norm_1 = norm(A2DRS_ATM_H, 1)

A2DRS_Basic_time = @elapsed begin
    A2DRS_Basic_H = a2drs_basic(A, lambda, problem, eps_opt)
end
A2DRS_Basic_H_norm_0 = matrix_norm_0(A2DRS_Basic_H)
A2DRS_Basic_H_norm_1 = norm(A2DRS_Basic_H, 1)

A2DRS_Basic_SG_time = @elapsed begin
    A2DRS_Basic_SG_H = a2drs_basic_sg(A, lambda, problem, eps_opt)
end
A2DRS_Basic_SG_H_norm_0 = matrix_norm_0(A2DRS_Basic_SG_H)
A2DRS_Basic_SG_H_norm_1 = norm(A2DRS_Basic_SG_H, 1)

# A2DRS_Halpern_VS_time = @elapsed begin
#     A2DRS_Halpern_VS_H = a2drs_halpern_vs(A, lambda, problem, eps_opt)
# end
# A2DRS_Halpern_VS_H_norm_0 = matrix_norm_0(A2DRS_Halpern_VS_H)
# A2DRS_Halpern_VS_H_norm_1 = norm(A2DRS_Halpern_VS_H, 1)

# A2DRS_Halpern_FS_time = @elapsed begin
#     A2DRS_Halpern_FS_H = a2drs_halpern_fs(A, lambda, problem, eps_opt)
# end
# A2DRS_Halpern_FS_H_norm_0 = matrix_norm_0(A2DRS_Halpern_FS_H)
# A2DRS_Halpern_FS_H_norm_1 = norm(A2DRS_Halpern_FS_H, 1)

# ADMM_time = @elapsed begin
#     ADMM_H = admm_p1(A)
# end
# ADMM_H_norm_0 = matrix_norm_0(ADMM_H)
# ADMM_H_norm_1 = norm(ADMM_H, 1)

rho = 3.0
epsilon = 10^(-5)
eps_abs = epsilon
eps_rel = epsilon
fixed_tol = false
eps_opt = epsilon
time_limit = 2*60*60
eps_opt = epsilon

ADMM_time = @elapsed begin
    ADMM_H = admm_p123(A, rho, eps_abs, eps_rel, fixed_tol, eps_opt, time_limit)
end
ADMM_H_norm_0 = matrix_norm_0(ADMM_H)
ADMM_H_norm_1 = norm(ADMM_H, 1)

if problem == "P123"
    println("PLS violation: $(norm(A' * A * ADMM_H - A'))")
    println("P123 violation: $(norm(ADMM_H * A * pinv(A) - ADMM_H))")
    if norm(A' * A * ADMM_H - A') + norm(ADMM_H * A * pinv(A) - ADMM_H) > 10^(-5)
        throw(ErrorException("Failed Infeasibility Check."))
    end
end

println("GRB time: $GRB_time")
println("GRB norm 1: $GRB_H_norm_1")

println("DRS time: $DRS_time")
println("DRS norm 1: $DRS_H_norm_1")

# println("A2DRS_TR time: $A2DRS_time")
# println("A2DRS_TR norm 1: $A2DRS_H_norm_1")

println("A2DRS_Boyd time: $A2DRS_Boyd_time")
println("A2DRS_Boyd norm 1: $A2DRS_Boyd_H_norm_1")

# println("A2DRS_ATM time: $A2DRS_ATM_time")
# println("A2DRS_ATM norm 1: $A2DRS_ATM_H_norm_1")

println("A2DRS_Basic time: $A2DRS_Basic_time")
println("A2DRS_Basic norm 1: $A2DRS_Basic_H_norm_1")

println("A2DRS_Basic_SG time: $A2DRS_Basic_SG_time")
println("A2DRS_Basic_SG norm 1: $A2DRS_Basic_SG_H_norm_1")

# println("A2DRS_Halpern_VS time: $A2DRS_Halpern_VS_time")
# println("A2DRS_Halpern_VS norm 1: $A2DRS_Halpern_VS_H_norm_1")

# println("A2DRS_Halpern_FS time: $A2DRS_Halpern_FS_time")
# println("A2DRS_Halpern_FS norm 1: $A2DRS_Halpern_FS_H_norm_1")

println("ADMM time: $ADMM_time")
println("ADMM norm 1: $ADMM_H_norm_1")