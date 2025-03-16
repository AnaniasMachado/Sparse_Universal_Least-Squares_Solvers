

function get_k0(F_matrix::Matrix{Float64})
    m, n = size(F_matrix)
    min_norm = Inf
    min_norm_idx = -1
    for j in 1:n
        current_norm = norm(F_matrix[:, j])
        if current_norm < min_norm
            min_norm = current_norm
            min_norm_idx = j
        elseif current_norm == min_norm
            min_norm_idx = j
        end
    end
    return min_norm_idx
end

# function compute_alpha(lambda::Float64, k0::Int64, res_f::Vector)
#     m, n = size(res_f[1])
#     p = size(res_f)[1] - 1

#     model = Model(Gurobi.Optimizer)

#     @variable(model, alpha[1:p])

#     idx_list = [i for i in 1:size(res_f)[1] if i != k0]
#     A = res_f[k0] + sum(alpha[i] * (res_f[idx_list[i]] - res_f[k0]) for i in 1:p)
#     @objective(model, Min, sum((A[i, j])^2 for i in 1:m for j in 1:n) + lambda * sum((alpha[i])^2 for i in 1:p))

#     set_optimizer_attribute(model, "LogToConsole", 0)

#     optimize!(model)

#     status = termination_status(model)
#     if status == MOI.OPTIMAL
#         alpha_star = [value(alpha[i]) for i in 1:p]
#         return alpha_star
#     else
#         println("Status: $status")
#         throw(ErrorException("Model was not optimized successfully."))
#     end
# end

function compute_alpha(t::Float64, k0::Int64, F_matrix::Matrix{Float64})
    m, n = size(F_matrix)
    p = n - 1
    idx_list = [j for j in 1:n if j != k0]
    B = zeros(m, p)
    for j in 1:p
        B[:, j] = F_matrix[:, idx_list[j]] - F_matrix[:, k0]
    end
    b = - B' * F_matrix[:, k0]
    A =  t * I(p) + B' * B
    alpha = A \ b
    return alpha
end

function compute_r(gamma::Float64, k0::Int64, F_matrix::Matrix{Float64})
    m, n = size(F_matrix)
    p = n - 1
    r = (1 - gamma * m) * norm(F_matrix[:, k0])
    for j in 1:n
        if j != k0
            r += gamma * norm(F_matrix[:, j])
        end
    end
    return r
end

function compute_ghat(alpha::Vector{Float64}, k0::Int64, V_matrix::Matrix{Float64})
    m, n = size(V_matrix)
    ghat = V_matrix[:, k0]
    idx_list = [i for i in 1:n if i != k0]
    p = n - 1
    for i in 1:p
        ghat += alpha[i] * (V_matrix[:, idx_list[i]] - V_matrix[:, k0])
    end
    return ghat
end

function compute_ared(lambda::Float64, fpi::Function, alpha::Vector{Float64}, gamma::Float64, k0::Int64, V_matrix::Matrix{Float64}, F_matrix::Matrix{Float64})
    r = compute_r(gamma, k0, F_matrix)
    ghat = compute_ghat(alpha, k0, V_matrix)
    fghat = fpi(ghat, lambda) - ghat
    return r - norm(fghat)
end

function compute_pred(alpha::Vector{Float64}, c::Float64, gamma::Float64, k0::Int64, F_matrix::Matrix{Float64})
    m, n = size(F_matrix)
    idx_list = [i for i in 1:n if i != k0]
    p = n - 1
    r = compute_r(gamma, k0, F_matrix)
    fhat = F_matrix[:, k0]
    for i in 1:p
        fhat += alpha[i] * (F_matrix[:, idx_list[i]] - F_matrix[:, k0])
    end
    return r - c * norm(fhat)
end

function compute_rho(lambda::Float64, fpi::Function, alpha::Vector{Float64}, c::Float64, gamma::Float64, k0::Int64, V_matrix::Matrix{Float64}, F_matrix::Matrix{Float64})
    ared = compute_ared(lambda, fpi, alpha, gamma, k0, V_matrix, F_matrix)
    pred = compute_pred(alpha, c, gamma, k0, F_matrix)
    return ared / (pred + 10^(-5))
end

function update_mu(mu::Float64, p1::Float64, p2::Float64, eta1::Float64, eta2::Float64, rho::Float64)
    if rho < p1
        return eta1 * mu
    elseif (rho >= p1) && (rho <= p2)
        return mu
    elseif rho > p2
        return eta2 * mu
    else
        println("Rho: $rho")
        println("p1: $p1")
        println("p2: $p2")
        throw(ErrorException("Invalid value of rho, p1 or p2."))
    end
end

function a2drs_tr(A::Matrix{Float64}, lambda::Float64, problem::String, eps_opt::Float64)
    # Initial data
    m, n = size(A)
    Xh = zeros(n, m)
    Vk = zeros(n, m)
    # V = rand(n, m)
    # V = generalized_inverse(A)
    # V = pinv(A)
    # Projection data
    # println("Point 0")
    proj_data = get_proj_data(A , problem)
    V = proj_data.AMP
    # println("Point 1")
    # Anderson acceleration data
    fpi = drs_tr_fpi(A, proj_data, problem)
    mu = 1.0
    # p1 = 0.01
    p1 = 10^(-2)
    p2 = 0.25
    # eta1 = 2.0
    eta1 = 1.05
    eta2 = 0.25
    gamma = 10^(-4)
    c = 1.0
    M = 5
    F_matrix = Matrix{Float64}(undef, m * n, 0)
    V_matrix = Matrix{Float64}(undef, m * n, 0)
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
        end
        V += X - Xh
        pri_res = primal_residual(A, Xh, proj_data, problem)
        dual_res = dual_residual(A, Xh, V, lambda, proj_data, problem)
        if (pri_res <= eps_opt) && (dual_res <= eps_opt)
            println("DRS TR Convergence: k=$k")
            break
        end
        println("Iteration TR k: $k")
        println("Primal TR residual: $pri_res")
        println("Dual TR residual: $dual_res")
        # Anderson acceleration
        f = vec(V - Vk)
        F_matrix = hcat(F_matrix, f)
        V_matrix = hcat(V_matrix, vec(V))
        if size(F_matrix)[2] > M
            F_matrix = F_matrix[:, 2:end]
            V_matrix = V_matrix[:, 2:end]
        end
        k0 = get_k0(F_matrix)
        # t = mu * (norm(F_matrix[:, k0])^2)
        t = (norm(F_matrix[:, k0])^2) * 10^(-3)
        alpha = compute_alpha(t, k0, F_matrix)
        rho = compute_rho(lambda, fpi, alpha, c, gamma, k0, V_matrix, F_matrix)
        # mu = update_mu(mu, p1, p2, eta1, eta2, rho)
        if rho > p1
            ghat = compute_ghat(alpha, k0, list_V)
            V = reshape(ghat, n, m)
        else
            V = reshape(fpi(V_matrix[:, k0], lambda), n, m)
        end
        # println("Current t: $t")
        # println("Current mu: $mu")
        # println("Current Fk0_F: $(norm(F_matrix[:, k0]))")
    end
    return Xh, k
end