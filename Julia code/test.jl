using CSV
using DataFrames

include("types.jl")
include("utility.jl")
include("solvers.jl")
# include("drs.jl")
include("admm.jl")
include("./src/call_admm.jl")

m = 10
n = 5
r = 2
A = gen_random_rank_r_matrix(m, n, r)
A_norm_0 = matrix_norm_0(A)
A_norm_1 = norm(A, 1)

AMP = pinv(A)
AMP_norm_0 = matrix_norm_0(AMP)
AMP_norm_1 = norm(AMP, 1)

data = DataInst(A, m, n, r)
constraints = ["P1"]
opt_tol = 10^(-8)

GRB_time = @elapsed begin
    GRB_H = gurobi_solver(data, constraints, opt_tol)
end
GRB_H_norm_0 = matrix_norm_0(GRB_H)
GRB_H_norm_1 = norm(GRB_H, 1)

iter = 5000000

# DRS_time = @elapsed begin
#     DRS_H = drs(A, 0.2)
# end
# DRS_H_norm_0 = matrix_norm_0(DRS_H)
# DRS_H_norm_1 = norm(DRS_H, 1)

# A2DR_time = @elapsed begin
#     A2DR_H = drs(A, 0.2)
# end
# A2DR_H_norm_0 = matrix_norm_0(A2DR_H)
# A2DR_H_norm_1 = norm(A2DR_H, 1)

# ADMM_time = @elapsed begin
#     ADMM_H = admm_p1(A)
# end
# ADMM_H_norm_0 = matrix_norm_0(ADMM_H)
# ADMM_H_norm_1 = norm(ADMM_H, 1)

rho = 4.0
epsilon = 10^(-8)
time_limit = 2*60*60

ADMM_time = @elapsed begin
    ADMM_H = admm_p123(A, rho, epsilon, epsilon, time_limit)
end
ADMM_H_norm_0 = matrix_norm_0(ADMM_H)
ADMM_H_norm_1 = norm(ADMM_H, 1)

eps_opt = 10^(-8)
fixed_tol = false

ADMM_GAB_time = @elapsed begin
    ADMM_GAB_H = run_admm_p123(A, m, n, r, rho, epsilon, epsilon, eps_opt, fixed_tol)
end
ADMM_GAB_H_norm_0 = matrix_norm_0(ADMM_GAB_H.H)
ADMM_GAB_H_norm_1 = norm(ADMM_GAB_H.H, 1)

println(norm(GRB_H, 1))
println("GRB time: $GRB_time")
# println(norm(DRS_H, 1))
# println("DRS time: $DRS_time")
# println(norm(A2DR_H, 1))
# println("A2DR time: $A2DR_time")
println(norm(ADMM_H, 1))
println("ADMM time: $ADMM_time")
println(norm(ADMM_GAB_H.H, 1))
println("ADMM time: $ADMM_GAB_time")

# results = DataFrame(
#     m = [m],
#     n = [n],
#     r = [r],
#     A_norm_0 = [A_norm_0],
#     A_norm_1 = [A_norm_1],
#     AMP_norm_0 = [AMP_norm_0],
#     AMP_norm_1 = [AMP_norm_1],
#     GRB_H_norm_0 = [GRB_H_norm_0],
#     GRB_H_norm_1 = [GRB_H_norm_1],
#     DRS_H_norm_0 = [DRS_H_norm_0],
#     DRS_H_norm_1 = [DRS_H_norm_1],
#     A2DR_H_norm_0 = [A2DR_H_norm_0],
#     A2DR_H_norm_1 = [A2DR_H_norm_1],
#     ADMM_H_norm_0 = [ADMM_H_norm_0],
#     ADMM_H_norm_1 = [ADMM_H_norm_1]
# )

# CSV.write("results.csv", results)