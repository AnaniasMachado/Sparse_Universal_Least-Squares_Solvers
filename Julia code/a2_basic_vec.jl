using LinearAlgebra

# function compute_gamma(g::Vector{Float64}, F_matrix::Matrix{Float64})
#     m, n = size(F_matrix)
#     eta = 10^(-8)

#     model = Model(Gurobi.Optimizer)

#     @variable(model, gamma[1:n])

#     FT = g - F_matrix * gamma
#     FT_F = sum((FT[i])^2 for i in 1:m)
#     ST = eta * sum((gamma[i])^2 for i in 1:n)
#     @objective(model, Min, FT_F + ST)

#     set_optimizer_attribute(model, "LogToConsole", 0)

#     optimize!(model)

#     status = termination_status(model)
#     if status == MOI.OPTIMAL
#         gamma_star = [value(gamma[i]) for i in 1:n]
#         return gamma_star
#     else
#         println("Status: $status")
#         throw(ErrorException("Model was not optimized successfully."))
#     end
# end

# function compute_initial_B(F_matrix::Matrix{Float64})
#     m, n = size(F_matrix)
#     eta = 10^(-2)
#     B = eta * I(n) + F_matrix' * F_matrix
#     return B
# end

# function update_B_F(kmod::Int64, F_matrix::Matrix{Float64}, f::Vector{Float64}, B::Matrix{Float64})
#     m, n = size(F_matrix)
#     u = f - F_matrix[:, kmod]
#     v = zeros(n)
#     v[kmod] = 1
#     uvT = u * v'
#     # B += F_matrix' * uvT + uvT' * F_matrix + u' * u * v * v'
#     B += F_matrix' * u * v' + v * u' * F_matrix + u' * u * v * v'
#     F_matrix[:, kmod] = f
#     return B, F_matrix
# end

function compute_gamma(g::Vector{Float64}, F_matrix::Matrix{Float64})
    m, n = size(F_matrix)
    eta = 10^(-2)

    A = eta * I(n) + F_matrix' * F_matrix
    b = F_matrix' * g
    return A \ b
end

# function compute_gamma(g::Vector{Float64}, F_matrix::Matrix{Float64}, B::Matrix{Float64})
#     b = F_matrix' * g
#     return B \ b
# end

# function compute_alpha(gamma::Vector{Float64})
#     p = size(gamma)[1]
#     alpha = zeros(p)
#     for i in 1:p
#         if i == 1
#             alpha[i] = gamma[i]
#         elseif (i > 1) && (i < p)
#             alpha[i] = gamma[i] - gamma[i-1]
#         else
#             alpha[i] = 1 - gamma[i-1]
#         end
#     end
#     return alpha
# end

function compute_alpha(gamma::Vector{Float64}, kmod::Int64)
    p = size(gamma)[1]
    alpha = zeros(p)
    idx_list = []
    for i in 0:p-1
        j = ((kmod + i) % p) + 1
        push!(idx_list, j)
    end
    for i in 1:p
        j = idx_list[i]
        if i == 1
            alpha[j] = gamma[j]
        elseif i == p
            jm = idx_list[i-1]
            alpha[j] = 1 - gamma[jm]
        else
            jm = idx_list[i-1]
            alpha[j] = gamma[j] - gamma[jm]
        end
    end
    return alpha
end

function compute_V_AA(V_matrix::Matrix{Float64}, alpha::Vector{Float64})
    return V_matrix * alpha
end

# function a2drs_basic(A::Matrix{Float64}, lambda::Float64, problem::String, eps_opt::Float64)
#     # Initial data
#     m, n = size(A)
#     Xh = rand(n, m)
#     # V = rand(n, m)
#     # V = pinv(A)
#     # V = generalized_inverse(A)
#     # Projection data
#     proj_data = get_proj_data(A , problem)
#     # Anderson acceleration data
#     fpi = drs_fpi(A, proj_data, problem)
#     M_max = 10
#     F_matrix = Matrix{Float64}(undef, m * n, 0)
#     V_matrix = Matrix{Float64}(undef, m * n, 0)
#     # Vk = rand(n, m)
#     # Vk = zeros(n, m)
#     # Vk = generalized_inverse(A)
#     Vk = proj_data.AMP
#     V = fpi(Vk, lambda)
#     gk = vec(Vk - V)
#     g = rand(n * m)
#     Vk = V
#     k = 0
#     while true
#         k += 1
#         # Usual DRS steps
#         Xh = soft_thresholding_matrix(V, lambda)
#         Vh = 2 * Xh - V
#         X = projection(A, Vh, proj_data, problem)
#         # X = gurobi_projection(Vh, proj_data, problem)
#         if !is_feasible(A, X, proj_data, problem)
#             println("Infeasible X.")
#             throw(ErrorException("Infeasible X error."))
#             break
#         else
#             # println("Feasible X.")
#         end
#         V += X - Xh
#         pri_res = primal_residual(A, Xh, proj_data, problem)
#         dual_res = dual_residual(A, Xh, V, lambda, proj_data, problem)
#         if (pri_res <= eps_opt) && (dual_res <= eps_opt)
#             println("DRS Basic Convergence: k=$k")
#             break
#         end
#         println("Iteration Basic k: $k")
#         println("Primal Basic residual: $pri_res")
#         println("Dual Basic residual: $dual_res")
#         # Anderson acceleration steps
#         ## Memory update
#         M = min(k, M_max)
#         V_DRS = vec(V)
#         g = vec(Vk) - V_DRS
#         f = g - gk
#         F_matrix = hcat(F_matrix, f)
#         V_matrix = hcat(V_matrix, V_DRS)
#         if size(F_matrix)[2] > M
#             F_matrix = F_matrix[:, 2:end]
#             V_matrix = V_matrix[:, 2:end]
#         end
#         if k == M_max
#             ## AA candidate
#             gamma = compute_gamma(g, F_matrix)
#             alpha = compute_alpha(gamma)
#             V_AA = compute_V_AA(V_matrix, alpha)
#             V_AA = reshape(V_AA, n, m)
#             ## Swaping
#             Vk = V
#             gk = g
#             V = V_AA
#         elseif M == M_max
#             ## AA candidate
#             gamma = compute_gamma(g, F_matrix)
#             alpha = compute_alpha(gamma)
#             V_AA = compute_V_AA(V_matrix, alpha)
#             V_AA = reshape(V_AA, n, m)
#             ## Swaping
#             Vk = V
#             gk = g
#             V = V_AA
#         else
#             ## Swaping
#             Vk = V
#             gk = g
#         end
#     end
#     return Xh, k
# end

function a2drs_basic(A::Matrix{Float64}, lambda::Float64, problem::String, eps_opt::Float64)
    # Initial data
    m, n = size(A)
    Xh = rand(n, m)
    # V = rand(n, m)
    # V = pinv(A)
    # V = generalized_inverse(A)
    # Projection data
    proj_data = get_proj_data(A , problem)
    # Anderson acceleration data
    fpi = drs_fpi(A, proj_data, problem)
    M_max = 10
    F_matrix = zeros(n*m, M_max)
    V_matrix = zeros(n*m, M_max)
    # Vk = rand(n, m)
    # Vk = zeros(n, m)
    # Vk = generalized_inverse(A)
    Vk = proj_data.AMP
    V = fpi(Vk, lambda)
    gk = vec(Vk - V)
    g = rand(n * m)
    Vk = V
    k = 0
    while true
        k += 1
        # Usual DRS steps
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
        pri_res = primal_residual(A, Xh, proj_data, problem)
        dual_res = dual_residual(A, Xh, V, lambda, proj_data, problem)
        if (pri_res <= eps_opt) && (dual_res <= eps_opt)
            println("DRS Basic Convergence: k=$k")
            break
        end
        println("Iteration Basic k: $k")
        println("Primal Basic residual: $pri_res")
        println("Dual Basic residual: $dual_res")
        # Anderson acceleration steps
        ## Memory update
        M = min(k, M_max)
        V_DRS = vec(V)
        g = vec(Vk) - V_DRS
        f = g - gk
        kmod = (k % M_max) + 1
        F_matrix[:, kmod] = f
        V_matrix[:, kmod] = V_DRS
        if k == M_max
            ## AA candidate
            gamma = compute_gamma(g, F_matrix)
            alpha = compute_alpha(gamma, kmod)
            V_AA = compute_V_AA(V_matrix, alpha)
            V_AA = reshape(V_AA, n, m)
            ## Swaping
            Vk = V
            gk = g
            V = V_AA
        elseif M == M_max
            ## AA candidate
            gamma = compute_gamma(g, F_matrix)
            alpha = compute_alpha(gamma, kmod)
            V_AA = compute_V_AA(V_matrix, alpha)
            V_AA = reshape(V_AA, n, m)
            ## Swaping
            Vk = V
            gk = g
            V = V_AA
        else
            ## Swaping
            Vk = V
            gk = g
        end
    end
    return Xh, k
end

function a2drs_basic_sg(A::Matrix{Float64}, lambda::Float64, problem::String, eps_opt::Float64)
    # Initial data
    m, n = size(A)
    Xh = rand(n, m)
    # V = rand(n, m)
    # V = pinv(A)
    # V = generalized_inverse(A)
    # Projection data
    proj_data = get_proj_data(A , problem)
    # Anderson acceleration data
    fpi = drs_fpi(A, proj_data, problem)
    M_max = 5
    F_matrix = Matrix{Float64}(undef, m * n, 0)
    V_matrix = Matrix{Float64}(undef, m * n, 0)
    # Vk = rand(n, m)
    # Vk = zeros(n, m)
    # Vk = generalized_inverse(A)
    Vkm = pinv(A)
    Vk = fpi(Vkm, lambda)
    V = rand(n, m)
    gkm = vec(Vkm - Vk)
    g0_F = norm(gkm)
    gk = rand(n * m)
    g = rand(n * m)
    Vk = V
    # Safeguard constants and variables
    D = 10^6
    R = 10
    n_AA = 0
    R_AA = 0
    I_sg = true
    k = 0
    while true
        k += 1
        # Usual DRS steps
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
        pri_res = primal_residual(A, Xh, proj_data, problem)
        dual_res = dual_residual(A, Xh, V, lambda, proj_data, problem)
        if (pri_res <= eps_opt) && (dual_res <= eps_opt)
            println("DRS Basic SG Convergence: k=$k")
            break
        end
        println("Iteration Basic SG k: $k")
        println("Primal Basic SG residual: $pri_res")
        println("Dual Basic SG residual: $dual_res")
        # Anderson acceleration steps
        ## Memory update
        M = min(k, M_max)
        V_DRS = vec(V)
        g = vec(Vk) - V_DRS
        f = vec(g - gk)
        F_matrix = hcat(F_matrix, f)
        V_matrix = hcat(V_matrix, V_DRS)
        if size(F_matrix)[2] > M
            F_matrix = F_matrix[:, 2:end]
            V_matrix = V_matrix[:, 2:end]
        end
        if k == M_max
            ## AA candidate
            gamma = compute_gamma(g, F_matrix)
            alpha = compute_alpha(gamma)
            V_AA = compute_V_AA(V_matrix, alpha)
            V_AA = reshape(V_AA, n, m)
            ## Swaping
            Vkm = Vk
            Vk = V
            gkm = gk
            gk = g
            V = V_AA
        elseif M == M_max
            ## AA candidate
            gamma = compute_gamma(g, F_matrix)
            alpha = compute_alpha(gamma)
            V_AA = compute_V_AA(V_matrix, alpha)
            V_AA = reshape(V_AA, n, m)
            ## Safeguard
            if I_sg || (R_AA >= R)
                if norm(gk) <= D * g0_F * (n_AA / (R + 1))^(-(1 + eps_opt))
                    V = V_AA
                    n_AA += 1
                    I_sg = false
                    R_AA = 1
                else
                    V = V_DRS
                    R_AA = 0
                end
            else
                V = V_AA
                n_AA += 1
                R_AA += 1
            end
            ## Swaping
            Vkm = Vk
            Vk = V
            gkm = gk
            gk = g
        else
            ## Swaping
            Vkm = Vk
            Vk = V
            gkm = gk
            gk = g
        end
    end
    return Xh, k
end