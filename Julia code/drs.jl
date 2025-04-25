using LinearAlgebra
using JuMP
using Gurobi

include("a2_tr.jl")

function soft_thresholding_matrix(X::Matrix{Float64}, lambda::Float64)
    return sign.(X) .* max.(abs.(X) .- lambda, 0)
end

function generalized_inverse(V::Matrix{Float64})
    m, n = size(V)

    model = Model(Gurobi.Optimizer)

    @variable(model, H[1:n, 1:m])

    @objective(model, Min, 0)

    @constraint(model, V * H * V .== V, base_name = "P1_")

    set_optimizer_attribute(model, "LogToConsole", 0)

    optimize!(model)

    status = termination_status(model)
    if status == MOI.OPTIMAL
        H_star = [value(H[i, j]) for i in 1:n, j in 1:m]
        return H_star
    else
        println("Model was not optimized successfully.")
    end
end

function gurobi_projection(V::Matrix{Float64}, proj_data::Union{DRSProjDataSimple, DRSProjDataP123}, problem::String)
    if problem == "P1"
        m, n = size(proj_data.T)

        model = Model(Gurobi.Optimizer)

        @variable(model, Y[1:n, 1:m])

        B = V - Y
        @objective(model, Min, sum((B[i, j])^2 for i in 1:size(B, 1) for j in 1:size(B, 2)))

        @constraint(model, proj_data.T * Y * proj_data.T .== proj_data.T, base_name = "P1_")

        set_optimizer_attribute(model, "LogToConsole", 0)

        optimize!(model)

        status = termination_status(model)
        if status == MOI.OPTIMAL
            Y_star = [value(Y[i, j]) for i in 1:n, j in 1:m]
            return Y_star
        else
            println("Model was not optimized successfully.")
        end
    end
end

# function get_proj_data(A::Matrix{Float64}, problem::String)
#     m, n = size(A)
#     AMP = pinv(A)
#     if problem == "P1"
#         R = A
#         S = A
#         T = A
#         println("Point 6")
#         RMTSM = AMP
#         SSMP = A * AMP
#         RMPR = AMP * A
#         println("Point 7")
#         V = 2 * (I(m * n) - kron(SSMP, RMPR))
#         println("Point 8")
#         VG = generalized_inverse(V)
#         println("Point 9")
#         proj_data = DRSProjData(R, S, T, RMTSM, SSMP, RMPR, VG)
#         println("Point 10")
#         return proj_data
#     end
# end

function get_proj_data(A::Matrix{Float64}, problem::String)
    m, n = size(A)
    if problem == "P1"
        AMP = pinv(A)
        R = A
        S = A
        T = A
        RMP = AMP
        SMP = AMP
        T_factor = RMP * T * SMP
        proj_data = DRSProjDataSimple(R, S, T, RMP, SMP, T_factor, AMP)
        return proj_data
    elseif problem == "PLS"
        R = A' * A
        S = I(m)
        T = A'
        RMP = pinv(R)
        SMP = I(m)
        T_factor = RMP * T
        AMP = RMP * A'
        proj_data = DRSProjDataSimple(R, S, T, RMP, SMP, T_factor, AMP)
        return proj_data
    elseif problem == "P123"
        AMP = pinv(A)
        proj_data = DRSProjDataP123(AMP)
        return proj_data
    end
end

# function projection(X::Matrix{Float64}, proj_data::DRSProjData)
#     Z = 2 * (proj_data.RMPR * X * proj_data.SSMP - X)
#     z = vec(Z)
#     w = proj_data.VG * z
#     W = reshape(w, size(X))
#     Y = proj_data.RMTSM + W - proj_data.RMPR * W * proj_data.SSMP
#     return Y
# end

function projection(A::Matrix{Float64}, X::Matrix{Float64}, proj_data::Union{DRSProjDataSimple, DRSProjDataP123}, problem::String)
    if problem == "P123"
        # L = 2 * proj_data.AMP * (proj_data.AMP' * X - proj_data.AMP' * proj_data.AMP)
        # Z = 2 * (proj_data.AMP * A * X - X)
        # G = Z - proj_data.AMP * A * Z
        # Y = X - 0.5 * (A' * A * L + proj_data.AMP * A * G - G)
        # return Y
        Z = proj_data.AMP * A * X * A * proj_data.AMP
        Y = proj_data.AMP - Z + X * A * proj_data.AMP
        return Y
    else
        Y = X - proj_data.RMP * proj_data.R * X * proj_data.S * proj_data.SMP + proj_data.T_factor
        return Y
    end
end

function primal_residual(A::Matrix, V::Matrix, proj_data::Union{DRSProjDataSimple, DRSProjDataP123}, problem::String)
    if problem == "P1"
        return norm(A * V * A - A)
    elseif problem == "PLS"
        return norm(A' * A * V - A')
    elseif problem == "P123"
        return norm(A' * V' * A' + V * A * proj_data.AMP - A' - V)
    end
end

function primal_residual_matrix(A::Matrix, V::Matrix, proj_data::Union{DRSProjDataSimple, DRSProjDataP123}, problem::String)
    if problem == "P1"
        return A * V * A - A
    elseif problem == "PLS"
        return A' * A * V - A'
    elseif problem == "P123"
        return A' * V' * A' + V * A * proj_data.AMP - A' - V
    end
end

function dual_variable(A::Matrix{Float64}, X::Matrix{Float64}, V::Matrix{Float64}, lambda::Float64, proj_data::Union{DRSProjDataSimple, DRSProjDataP123}, problem::String)
    if problem == "P123"
        B = (1 / lambda) * (X - V)
        L = proj_data.AMP * proj_data.AMP' * B * A * proj_data.AMP
        G = B * A * proj_data.AMP - B
        return L, G
    else
        L = proj_data.RMP' * (1 / lambda) * (X - V) * proj_data.SMP'
        return L
    end
end

function dual_residual(A::Matrix{Float64}, X::Matrix{Float64}, V::Matrix{Float64}, lambda::Float64, proj_data::Union{DRSProjDataSimple, DRSProjDataP123}, problem::String)
    if problem == "P123"
        Lambda, Gamma = dual_variable(A, X, V, lambda, proj_data, problem)
        return norm((1 / lambda) * (V - X) + A' * A * Lambda + Gamma * A * proj_data.AMP - Gamma)
    else
        Lambda = dual_variable(A, X, V, lambda, proj_data, problem)
        return norm((1 / lambda) * (V - X) + proj_data.R' * Lambda * proj_data.S')
    end
end

function dual_residual_matrix(A::Matrix{Float64}, X::Matrix{Float64}, V::Matrix{Float64}, lambda::Float64, proj_data::Union{DRSProjDataSimple, DRSProjDataP123}, problem::String)
    if problem == "P123"
        Lambda, Gamma = dual_variable(A, X, V, lambda, proj_data, problem)
        return (1 / lambda) * (V - X) + A' * A * Lambda + Gamma * A * proj_data.AMP - Gamma
    else
        Lambda = dual_variable(A, X, V, lambda, proj_data, problem)
        return (1 / lambda) * (V - X) + proj_data.R' * Lambda * proj_data.S'
    end
end

function is_feasible(A::Matrix{Float64}, X::Matrix{Float64}, proj_data::Union{DRSProjDataSimple, DRSProjDataP123}, problem::String)
    if problem == "P1"
        return norm(A * X * A - A) < 10^(-5)
    elseif problem == "PLS"
        return norm(A' * A * X - A') < 10^(-5)
    elseif problem == "P123"
        # println("PLS violation: $(norm(A' * A * X - A'))")
        # println("P123 violation: $(norm(X * A * proj_data.AMP - X))")
        return (norm(A' * A * X - A') < 10^(-5)) && (norm(X * A * proj_data.AMP - X) < 10^(-5))
    end
end

# function drs(A::Matrix{Float64}, lambda::Float64, problem::String, eps_opt::Float64)
#     # Initial data
#     m, n = size(A)
#     Xh = zeros(n, m)
#     # V = rand(n, m)
#     # V = pinv(A)
#     # V = generalized_inverse(A)
#     # Projection data
#     proj_data = get_proj_data(A , problem)
#     c_F = 0.0
#     if problem == "P123"
#         c_F = norm(A')
#     end
#     V = proj_data.AMP
#     k = 0
#     while true
#         k += 1
#         Xh = soft_thresholding_matrix(V, lambda)
#         Vh = 2 * Xh - V
#         X = projection(A, Vh, proj_data, problem)
#         # X = gurobi_projection(Vh, proj_data, problem)
#         # if !is_feasible(A, X, proj_data, problem)
#         #     println("Infeasible X.")
#         #     throw(ErrorException("Infeasible X error."))
#         #     break
#         # else
#         #     # println("Feasible X.")
#         # end
#         V += X - Xh
#         pri_res = primal_residual(A, Xh, proj_data, problem)
#         dual_res = dual_residual(A, Xh, V, lambda, proj_data, problem)
#         if (pri_res <= eps_opt) && (dual_res <= eps_opt)
#             println("DRS Convergence: k=$k")
#             break
#         end
#         # println("Iteration k: $k")
#         # println("Primal residual: $pri_res")
#         # println("Dual residual: $dual_res")
#     end
#     return Xh, k
# end

function drs(A::Matrix{Float64}, lambda::Float64, eps_abs::Float64, eps_rel::Float64, problem::String, fixed_tol::Bool, eps_opt::Float64, stop_crit::String, time_limit::Int64)
    # Initial data
    m, n = size(A)
    Xh = zeros(n, m)
    X = zeros(n, m)
    # V = rand(n, m)
    # V = pinv(A)
    # V = generalized_inverse(A)
    # Projection data
    proj_data = get_proj_data(A , problem)
    V = proj_data.AMP
    initial_pri_res = primal_residual_matrix(A, Xh, proj_data, problem)
    initial_dual_res = dual_res = dual_residual_matrix(A, Xh, V, lambda, proj_data, problem)
    r0 = norm(hcat(initial_pri_res, initial_dual_res))
    eps_tol = eps_abs + eps_rel * r0
    start_time = time()
    k = 0
    while true
        k += 1
        Xh = soft_thresholding_matrix(V, lambda)
        Vh = 2 * Xh - V
        X = projection(A, Vh, proj_data, problem)
        # X = gurobi_projection(Vh, proj_data, problem)
        # if !is_feasible(A, X, proj_data, problem)
        #     println("Infeasible X.")
        #     throw(ErrorException("Infeasible X error."))
        #     break
        # else
        #     # println("Feasible X.")
        # end
        V += X - Xh
        if fixed_tol
            pri_res = primal_residual(A, Xh, proj_data, problem)
            dual_res = dual_residual(A, Xh, V, lambda, proj_data, problem)
            if (pri_res <= eps_opt) && (dual_res <= eps_opt)
                # println("DRS Convergence: k=$k")
                # println("rp: $pri_res")
                # println("rd: $dual_res")
                break
            end
            # println("Iteration k: $k")
            # println("Primal residual: $pri_res")
            # println("Dual residual: $dual_res")
        elseif !fixed_tol && (stop_crit == "Boyd")
            pri_res = primal_residual_matrix(A, Xh, proj_data, problem)
            dual_res = dual_residual_matrix(A, Xh, V, lambda, proj_data, problem)
            res = hcat(pri_res, dual_res)
            if norm(res) <= eps_tol
                # println("DRS Convergence: k=$k")
                # println("rp: $(norm(pri_res))")
                # println("rd: $(norm(dual_res))")
                break
            end
        elseif !fixed_tol && (stop_crit == "Fixed_Point")
            res = norm(X - Xh)
            if norm(res) <= eps_tol
                # println("DRS Convergence: k=$k")
                # println("res: $res")
                break
            end
        end
        # Checks time limit
        elapsed_time = time() - start_time
        if elapsed_time > time_limit
            println("TimeLimit: DRS exceed time limit to solve the problem.")
            return "-", k
        end
    end
    return X, k
end

function drs_fpi(A::Matrix{Float64}, proj_data::Union{DRSProjDataSimple, DRSProjDataP123}, problem::String)
    function fpi(V::Matrix{Float64}, lambda::Float64)
        Xh = soft_thresholding_matrix(V, lambda)
        Vh = 2 * Xh - V
        X = projection(A, Vh, proj_data, problem)
        return V + X - Xh
    end
    return fpi
end

function drs_tr_fpi(A::Matrix{Float64}, proj_data::Union{DRSProjDataSimple, DRSProjDataP123}, problem::String)
    m, n = size(A)
    function tr_fpi(V::Vector{Float64}, lambda::Float64)
        V = reshape(V, n, m)
        Xh = soft_thresholding_matrix(V, lambda)
        Vh = 2 * Xh - V
        X = projection(A, Vh, proj_data, problem)
        V += X - Xh
        return vec(V)
    end
    return tr_fpi
end