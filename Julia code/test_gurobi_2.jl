using MAT
using Base.GC

include("types.jl")
#include("utility.jl")
include("solvers_2.jl")

mat_filepath = "experiment_6_matrix_m180_n90_r30_d25_idx5.mat"
mat_data = matread(mat_filepath)
A = mat_data["matrix"]
A = Matrix(A)

m_value = match(r"m(\d+)", mat_filepath).captures[1]
n_value = match(r"n(\d+)", mat_filepath).captures[1]
r_value = match(r"r(\d+)", mat_filepath).captures[1]
d_value = match(r"d(\d+)", mat_filepath).captures[1]

m = parse(Int, m_value)
n = parse(Int, n_value)
r = parse(Int, r_value)
d = parse(Int, d_value)

hatA = A' * A
hatm, hatn = size(hatA)
data = DataInst(hatA, hatm, hatn, r)

constraints = ["P1"]
opt_tol = 10^(-5)

GRB_P1_time = @elapsed begin
    GRB_P1_H = gurobi_solver(data, constraints, opt_tol)
end
#GRB_P1_H_norm_0 = matrix_norm_0(GRB_P1_H)
#GRB_P1_H_norm_1 = norm(GRB_P1_H, 1)

# GC.gc()

# constraints = ["P1", "Sym"]

# GRB_P1_Sym_time = @elapsed begin
#     GRB_P1_Sym_H = gurobi_solver(data, constraints, opt_tol)
# end
# #GRB_P1_Sym_H_norm_0 = matrix_norm_0(GRB_P1_Sym_H)
# #GRB_P1_Sym_H_norm_1 = norm(GRB_P1_Sym_H, 1)

# @show GRB_P1_time, GRB_P1_Sym_time
@show GRB_P1_time