

# function compute_gamma(g::Matrix{Float64}, F_list::Vector)
#     m, n = size(g)
#     p = size(F_list)[1]
#     eta = 0.2

#     model = Model(Gurobi.Optimizer)

#     @variable(model, gamma[1:p])

#     FT = g - sum(F_list[i] * gamma[i] for i in 1:p)
#     # FT_F = sum((FT[i, j])^2 for i in 1:m for j in 1:n)
#     # ST = sum((F_list[k][i, j])^2 for i in 1:m for j in 1:n for k in 1:p)
#     ST = eta * sum((gamma[i])^2 for i in 1:p)
#     # TT = eta * ST * sum((gamma[i])^2 for i in 1:p)
#     # @objective(model, Min, sum((FT[i, j])^2 for i in 1:m for j in 1:n))
#     @objective(model, Min, ST + sum((FT[i, j])^2 for i in 1:m for j in 1:n))
#     # @objective(model, Min, FT_F + TT)

#     set_optimizer_attribute(model, "LogToConsole", 0)

#     optimize!(model)

#     status = termination_status(model)
#     if status == MOI.OPTIMAL
#         gamma_star = [value(gamma[i]) for i in 1:p]
#         return gamma_star
#     else
#         println("Status: $status")
#         throw(ErrorException("Model was not optimized successfully."))
#     end
# end

# function compute_gamma(g::Matrix{Float64}, F_list::Vector)
#     m, n = size(g)
#     p = size(F_list)[1]
#     eta = 0.2
#     B = zeros(p, p)
#     b = zeros(p)
#     for i in 1:p
#         b[i] = tr(g' * F_list[i])
#         for j in i:p
#             t = tr(F_list[j]' * F_list[i])
#             B[i, j] = t
#             B[j, i] = t
#         end
#     end
#     A = 2 * (eta * I(p) + B)
#     gamma = A \ b
#     return gamma
# end

function compute_gamma(B::Matrix{Float64}, b::Vector{Float64})
    p = size(b)[1]
    eta = 0.2
    A = 2 * (eta * I(p) + B)
    gamma = A \ b
    return gamma
end

function compute_alpha(gamma::Vector{Float64})
    p = size(gamma)[1]
    alpha = zeros(p)
    for i in 1:p
        if i == 1
            alpha[i] = gamma[i]
        elseif (i > 1) && (i < p)
            alpha[i] = gamma[i] - gamma[i-1]
        else
            alpha[i] = 1 - gamma[i-1]
        end
    end
    return alpha
end

function compute_V_AA(V_list::Vector, alpha::Vector{Float64})
    m, n = size(V_list[1])
    p = size(alpha)[1]
    V_AA = zeros(m, n)
    for i in 1:p
        V_AA += alpha[i] * V_list[i]
    end
    return V_AA
end

function compute_initial_B_b(g::Matrix{Float64}, F_list::Vector)
    p = size(F_list)[1]
    B = zeros(p, p)
    b = zeros(p)
    for i in 1:p
        b[i] = tr(g' * F_list[i])
        for j in i:p
            t = tr(F_list[j]' * F_list[i])
            B[i, j] = t
            B[j, i] = t
        end
    end
    return B, b
end

function update_B_b(B::Matrix{Float64}, b::Vector{Float64}, g::Matrix{Float64}, F_list::Vector)
    p = size(F_list)[1]
    B[1:p-1, 1:p-1] = B[2:p, 2:p]
    b[1:p-1] = b[2:p]
    for i in 1:p
        t = tr(F_list[i]' * F_list[p])
        B[i, p] = t
        B[p, i] = t
    end
    b[p] = tr(g' * F_list[p])
    return B, b
end

function a2drs_basic(A::Matrix{Float64}, lambda::Float64, problem::String, eps_opt::Float64)
    # Initial data
    m, n = size(A)
    X = rand(n, m)
    # V = rand(n, m)
    # V = pinv(A)
    # V = generalized_inverse(A)
    # Projection data
    proj_data = get_proj_data(A , problem)
    # Anderson acceleration data
    fpi = drs_fpi(A, proj_data, problem)
    M_max = 5
    F_list = []
    V_list = []
    # Vk = rand(n, m)
    # Vk = zeros(n, m)
    # Vk = generalized_inverse(A)
    Vk = pinv(A)
    V = fpi(Vk, lambda)
    gk = Vk - V
    g = rand(n, m)
    Vk = V
    B = zeros(M_max, M_max)
    b = zeros(M_max)
    k = 0
    while true
        k += 1
        # Usual DRS steps
        Xh = soft_thresholding_matrix(V, lambda)
        Vh = 2 * Xh - V
        X = projection(A, Vh, proj_data, problem)
        # X = gurobi_projection(Vh, proj_data, problem)
        if !is_feasible(A, X, proj_data, problem)
            println("Infeasible X.")
            throw(ErrorException("Infeasible X error."))
            break
        else
            # println("Feasible X.")
        end
        V += X - Xh
        pri_res = primal_residual(A, Xh, proj_data, problem)
        dual_res = dual_residual(A, Xh, V, lambda, proj_data, problem)
        if (pri_res <= eps_opt) && (dual_res <= eps_opt)
            println("DRS Basic Convergence: k=$k")
            break
        end
        # println("Iteration Basic k: $k")
        # println("Primal Basic residual: $pri_res")
        # println("Dual Basic residual: $dual_res")
        # Anderson acceleration steps
        ## Memory update
        M = min(k, M_max)
        V_DRS = V
        g = Vk - V_DRS
        f = g - gk
        push!(F_list, f)
        push!(V_list, V_DRS)
        if size(F_list)[1] > M
            splice!(F_list, 1)
            splice!(V_list, 1)
        end
        if k == M_max
            B, b = compute_initial_B_b(g, F_list)
            ## AA candidate
            gamma = compute_gamma(B, b)
            alpha = compute_alpha(gamma)
            V_AA = compute_V_AA(V_list, alpha)
            ## Swaping
            Vk = V
            gk = g
            V = V_AA
        elseif M == M_max
            B, b = update_B_b(B, b, g, F_list)
            ## AA candidate
            gamma = compute_gamma(B, b)
            alpha = compute_alpha(gamma)
            V_AA = compute_V_AA(V_list, alpha)
            ## Swaping
            Vk = V
            gk = g
            V = V_AA
        else
            ## Swaping
            Vk = V
            gk = g
        end
        # g_AA = Vk - V_AA
        # if norm(g_AA) < norm(g)
        #     V = V_AA
        # else
        #     V = V_DRS
        # end
    end
    return X
end

# function a2drs_basic_sg(A::Matrix{Float64}, lambda::Float64, problem::String, eps_opt::Float64)
#     # Initial data
#     m, n = size(A)
#     X = rand(n, m)
#     # V = rand(n, m)
#     # V = pinv(A)
#     # V = generalized_inverse(A)
#     # Projection data
#     proj_data = get_proj_data(A , problem)
#     # Anderson acceleration data
#     fpi = drs_fpi(A, proj_data, problem)
#     M_max = 5
#     F_list = []
#     V_list = []
#     # Vk = rand(n, m)
#     # Vk = zeros(n, m)
#     # Vk = generalized_inverse(A)
#     Vk = pinv(A)
#     V = fpi(Vk, lambda)
#     gk = Vk - V
#     g0_F = norm(gk)
#     g = rand(n, m)
#     Vk = V
#     B = zeros(M_max, M_max)
#     b = zeros(M_max)
#     # Safeguard constants and variables
#     D = 10^6
#     R = 10
#     n_AA = 0
#     R_AA = 0
#     I_sg = true
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
#             println("DRS Basic SG Convergence: k=$k")
#             break
#         end
#         # println("Iteration Basic k: $k")
#         # println("Primal Basic residual: $pri_res")
#         # println("Dual Basic residual: $dual_res")
#         # Anderson acceleration steps
#         ## Memory update
#         M = min(k, M_max)
#         V_DRS = V
#         g = Vk - V_DRS
#         f = g - gk
#         push!(F_list, f)
#         push!(V_list, V_DRS)
#         if size(F_list)[1] > M
#             splice!(F_list, 1)
#             splice!(V_list, 1)
#         end
#         if k == M_max
#             B, b = compute_initial_B_b(g, F_list)
#             ## AA candidate
#             gamma = compute_gamma(B, b)
#             alpha = compute_alpha(gamma)
#             V_AA = compute_V_AA(V_list, alpha)
#             ## Swaping
#             Vk = V
#             gk = g
#             V = V_AA
#         elseif M == M_max
#             B, b = update_B_b(B, b, g, F_list)
#             ## AA candidate
#             gamma = compute_gamma(B, b)
#             alpha = compute_alpha(gamma)
#             V_AA = compute_V_AA(V_list, alpha)
#             ## Safeguard
#             if I_sg || (R_AA >= R)
#                 if norm(gk) <= D * g0_F * (n_AA / (R + 1))^(-(1 + eps_opt))
#                     V = V_AA
#                     n_AA += 1
#                     I_sg = false
#                     R_AA = 1
#                 else
#                     V = V_DRS
#                     R_AA = 0
#                 end
#             else
#                 V = V_AA
#                 n_AA += 1
#                 R_AA += 1
#             end
#             ## Swaping
#             Vk = V
#             gk = g
#         else
#             ## Swaping
#             Vk = V
#             gk = g
#         end
#         # g_AA = Vk - V_AA
#         # if norm(g_AA) < norm(g)
#         #     V = V_AA
#         # else
#         #     V = V_DRS
#         # end
#     end
#     return X
# end

function a2drs_basic_sg(A::Matrix{Float64}, lambda::Float64, problem::String, eps_opt::Float64)
    # Initial data
    m, n = size(A)
    X = rand(n, m)
    # V = rand(n, m)
    # V = pinv(A)
    # V = generalized_inverse(A)
    # Projection data
    proj_data = get_proj_data(A , problem)
    # Anderson acceleration data
    fpi = drs_fpi(A, proj_data, problem)
    M_max = 5
    F_list = []
    V_list = []
    # Vk = rand(n, m)
    # Vk = zeros(n, m)
    # Vk = generalized_inverse(A)
    Vkm = pinv(A)
    Vk = fpi(Vkm, lambda)
    V = rand(n, m)
    gkm = Vkm - Vk
    g0_F = norm(gkm)
    gk = rand(n, m)
    g = rand(n, m)
    Vk = V
    B = zeros(M_max, M_max)
    b = zeros(M_max)
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
        if !is_feasible(A, X, proj_data, problem)
            println("Infeasible X.")
            throw(ErrorException("Infeasible X error."))
            break
        else
            # println("Feasible X.")
        end
        V += X - Xh
        pri_res = primal_residual(A, Xh, proj_data, problem)
        dual_res = dual_residual(A, Xh, V, lambda, proj_data, problem)
        if (pri_res <= eps_opt) && (dual_res <= eps_opt)
            println("DRS Basic SG Convergence: k=$k")
            break
        end
        # println("Iteration Basic k: $k")
        # println("Primal Basic residual: $pri_res")
        # println("Dual Basic residual: $dual_res")
        # Anderson acceleration steps
        ## Memory update
        M = min(k, M_max)
        V_DRS = V
        g = Vk - V_DRS
        f = gk - gkm
        push!(F_list, f)
        push!(V_list, V_DRS)
        if size(F_list)[1] > M
            splice!(F_list, 1)
            splice!(V_list, 1)
        end
        if k == M_max
            B, b = compute_initial_B_b(g, F_list)
            ## AA candidate
            gamma = compute_gamma(B, b)
            alpha = compute_alpha(gamma)
            V_AA = compute_V_AA(V_list, alpha)
            ## Swaping
            Vkm = Vk
            Vk = V
            gkm = gk
            gk = g
            V = V_AA
        elseif M == M_max
            B, b = update_B_b(B, b, g, F_list)
            ## AA candidate
            gamma = compute_gamma(B, b)
            alpha = compute_alpha(gamma)
            V_AA = compute_V_AA(V_list, alpha)
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
        # g_AA = Vk - V_AA
        # if norm(g_AA) < norm(g)
        #     V = V_AA
        # else
        #     V = V_DRS
        # end
    end
    return X
end