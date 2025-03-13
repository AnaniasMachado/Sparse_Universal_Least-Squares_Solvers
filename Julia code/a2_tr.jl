

function get_k0(res_f::Vector)
    min_norm = Inf
    min_norm_idx = -1
    for (i, f) in enumerate(res_f)
        current_norm = norm(f)
        if current_norm < min_norm
            min_norm = current_norm
            min_norm_idx = i
        elseif current_norm == min_norm
            min_norm_idx = max(min_norm_idx, i)
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

function compute_alpha(lambda::Float64, k0::Int64, res_f::Vector)
    m, n = size(res_f[1])
    p = size(res_f)[1] - 1
    B = zeros(p, p)
    b = zeros(p)
    idx_list = [i for i in 1:size(res_f)[1] if i != k0]

    for i in 1:p
        b[i] = -tr(res_f[k0]' * (res_f[idx_list[i]] - res_f[k0]))
        for j in 1:p
            t = tr((res_f[idx_list[j]] - res_f[k0])' * (res_f[idx_list[i]] - res_f[k0]))
            B[i, j] = t
            B[j, i] = t
        end
    end

    A = 2 * (lambda * I(p) + B)
    alpha = A \ b
    return alpha
end

function compute_r(gamma::Float64, k0::Int64, res_f::Vector)
    m = size(res_f)[1] - 1
    r = (1 - gamma * m) * norm(res_f[k0])
    for i in 1:size(res_f)[1]
        if i != k0
            r += gamma * norm(res_f[i])
        end
    end
    return r
end

function compute_ghat(alpha::Vector{Float64}, k0::Int64, list_V::Vector)
    ghat = list_V[k0]
    idx_list = [i for i in 1:size(list_V)[1] if i != k0]
    p = size(list_V)[1] - 1
    for i in 1:p
        ghat += alpha[i] * (list_V[idx_list[i]] - list_V[k0])
    end
    return ghat
end

# function compute_ared(t::Float64, alpha::Vector{Float64}, gamma::Float64, k0::Int64, list_V::Vector, res_f::Vector)
function compute_ared(t::Float64, fpi::Function, alpha::Vector{Float64}, gamma::Float64, k0::Int64, list_V::Vector, res_f::Vector)
    r = compute_r(gamma, k0, res_f)
    ghat = compute_ghat(alpha, k0, list_V)
    # fghat = soft_thresholding_matrix(ghat, t) - ghat
    fghat = fpi(ghat, t) - ghat
    return r - norm(fghat)
end

function compute_pred(alpha::Vector{Float64}, c::Float64, gamma::Float64, k0::Int64, res_f::Vector)
    r = compute_r(gamma, k0, res_f)
    fhat = res_f[k0]
    idx_list = [i for i in 1:size(res_f)[1] if i != k0]
    p = size(res_f)[1] - 1
    for i in 1:p
        fhat += alpha[i] * (res_f[idx_list[i]] - res_f[k0])
    end
    return r - c * norm(fhat)
end

# function compute_rho(t::Float64, alpha::Vector{Float64}, c::Float64, gamma::Float64, k0::Int64, list_V::Vector, res_f::Vector)
function compute_rho(t::Float64, fpi::Function, alpha::Vector{Float64}, c::Float64, gamma::Float64, k0::Int64, list_V::Vector, res_f::Vector)
    # ared = compute_ared(t, alpha, gamma, k0, list_V, res_f)
    ared = compute_ared(t, fpi, alpha, gamma, k0, list_V, res_f)
    pred = compute_pred(alpha, c, gamma, k0, res_f)
    # println("Current ared: $ared")
    # println("Current pred: $pred")
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

# function anderson_acceleration(V::Matrix{Float64}, t::Float64, mu::Float64, p1::Float64, p2::Float64, eta1::Float64, eta2::Float64, gamma::Float64, c::Float64, M::Int64, list_V::Vector, res_f::Vector, eps_opt::Float64)
function anderson_acceleration_tr(V::Matrix{Float64}, t::Float64, fpi::Function, mu::Float64, p1::Float64, p2::Float64, eta1::Float64, eta2::Float64, gamma::Float64, c::Float64, M::Int64, list_V::Vector, res_f::Vector, eps_opt::Float64)
    # g = soft_thresholding_matrix(V, t)
    g = fpi(V, t)
    f = g - V
    if norm(f) < eps_opt
        return g
    end
    push!(list_V, g)
    push!(res_f, f)
    if size(res_f)[1] == 1
        return g
    end
    if size(res_f)[1] > M
        splice!(list_V, 1)
        splice!(res_f, 1)
    end
    k0 = get_k0(res_f)
    lambda = mu * (norm(res_f[k0])^2)
    alpha = compute_alpha(lambda, k0, res_f)
    rho = compute_rho(t, fpi, alpha, c, gamma, k0, list_V, res_f)
    mu = update_mu(mu, p1, p2, eta1, eta2, rho)
    if rho > p1
        ghat = compute_ghat(alpha, k0, list_V)
        return ghat
    else
        return list_V[k0]
    end
end