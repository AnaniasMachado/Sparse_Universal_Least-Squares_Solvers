# function compute_gamma(g::Vector{Float64}, F_matrix::Matrix{Float64}, S_matrix::Matrix{Float64})
#     m, n = size(F_matrix)
#     eta = 10^(-8)

#     model = Model(Gurobi.Optimizer)

#     @variable(model, gamma[1:n])

#     FT = g - F_matrix * gamma
#     FT_F = sum((FT[i])^2 for i in 1:m)
#     ST = eta * (norm(F_matrix)^2 + norm(S_matrix)^2) * sum((gamma[i])^2 for i in 1:n)
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

function compute_gamma(g::Vector{Float64}, F_matrix::Matrix{Float64}, S_matrix::Matrix{Float64})
    m, n = size(F_matrix)
    eta = 10^(-8)

    A = eta * (norm(F_matrix)^2 + norm(S_matrix)^2) * I(n) + F_matrix' * F_matrix
    b = F_matrix' * g
    return A \ b
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

function compute_V_AA(V_matrix::Matrix{Float64}, alpha::Vector{Float64})
    return V_matrix * alpha
end

function a2drs_boyd(A::Matrix{Float64}, lambda::Float64, problem::String, eps_opt::Float64)
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
    F_matrix = Matrix{Float64}(undef, m * n, 0)
    V_matrix = Matrix{Float64}(undef, m * n, 0)
    S_matrix = Matrix{Float64}(undef, m * n, 0)
    # Vk = rand(n, m)
    # Vk = zeros(n, m)
    # Vk = generalized_inverse(A)
    Vkm = proj_data.AMP
    Vk = fpi(Vkm, lambda)
    Vkp = zeros(n, m)
    gkm = vec(Vkm - Vk)
    g0_F = norm(gkm)
    gk = rand(n * m)
    g = rand(n * m)
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
        Xh = soft_thresholding_matrix(Vk, lambda)
        Vh = 2 * Xh - Vk
        X = projection(A, Vh, proj_data, problem)
        # X = gurobi_projection(Vh, proj_data, problem)
        # if !is_feasible(A, X, proj_data, problem)
        #     println("Infeasible X.")
        #     throw(ErrorException("Infeasible X error."))
        #     break
        # else
        #     # println("Feasible X.")
        # end
        Vkp = Vk + X - Xh
        pri_res = primal_residual(A, Xh, proj_data, problem)
        dual_res = dual_residual(A, Xh, Vkp, lambda, proj_data, problem)
        if (pri_res <= eps_opt) && (dual_res <= eps_opt)
            println("A2DRS Boyd Convergence: k=$k")
            break
        end
        println("Iteration Boyd k: $k")
        println("Primal Boyd residual: $pri_res")
        println("Dual Boyd residual: $dual_res")
        # Anderson acceleration steps
        ## Memory update
        M = min(k, M_max)
        V_DRS = vec(Vkp)
        g = vec(Vk) - V_DRS
        f = vec(g - gk)
        s = vec(Vk - Vkm)
        F_matrix = hcat(F_matrix, f)
        V_matrix = hcat(V_matrix, V_DRS)
        S_matrix = hcat(S_matrix, s)
        if size(F_matrix)[2] > M
            F_matrix = F_matrix[:, 2:end]
            V_matrix = V_matrix[:, 2:end]
            S_matrix = S_matrix[:, 2:end]
        end
        ## AA candidate
        gamma = compute_gamma(g, F_matrix)
        alpha = compute_alpha(gamma)
        V_AA = compute_V_AA(V_matrix, alpha)
        V_AA = reshape(V_AA, n, m)
        ## Safeguard
        if I_sg || (R_AA >= R)
            if norm(gk) <= D * g0_F * (n_AA / (R + 1))^(-(1 + eps_opt))
                Vkp = V_AA
                n_AA += 1
                I_sg = false
                R_AA = 1
            else
                Vkp = V_DRS
                R_AA = 0
            end
        else
            Vkp = V_AA
            n_AA += 1
            R_AA += 1
        end
        ## Swaping
        Vkm = Vk
        Vk = Vkp
        gkm = gk
        gk = g
    end
    return Xh, k
end