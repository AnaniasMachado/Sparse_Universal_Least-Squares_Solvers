using LinearAlgebra

epsilon = 10^(-5)

function soft_thresholding_matrix(X::Matrix{Float64}, lambda::Float64)
    return sign.(X) .* max.(abs.(X) .- lambda, 0)
end

function count_singular_values(S::Diagonal)
    rank = count(x -> abs(x) > epsilon, S)
    return rank
end

function variables_initialization(V1::Matrix{Float64}, U1::Matrix{Float64}, D_inv::Matrix{Float64}, rho::Float64)
    Theta = (V1 * U1') / norm(V1 * U1', Inf)
    Lambda = Theta / rho
    E = V1 * D_inv * U1' + Lambda
    return Lambda, E
end

function admm_p123(A::Matrix{Float64}, rho::Float64, eps_abs::Float64, eps_rel::Float64, fixed_tol::Float64, time_limit::Int64)
    m, n = size(A)
    U, S, V = svd(A)
    S = Diagonal(S)
    r = count_singular_values(S)
    U1 = U[:, 1:r]
    V1 = V[:, 1:r]
    V2 = V[:, r+1:end]
    D_inv = zeros(r, r)
    for i in 1:r
        D_inv[i, i] = 1 / S[i, i]
    end
    Z = zeros(n - r, r)
    Lambda, Ekm = variables_initialization(V1, U1, D_inv, rho)
    start_time = time()
    while true
        # Updates J
        J = Ekm - V1 * D_inv * U1' - Lambda
        # Updates Zk
        Z = V2' * J * U1
        # Updates Y
        Y = V1 * D_inv * U1' + V2 * Z * U1' + Lambda
        # Updates Ek
        Ek = soft_thresholding_matrix(Y, 1/rho)
        # Updates Lambdak
        Lambda = Lambda + V1 * D_inv * U1' + V2 * Z * U1' - Ek
        # Calculates stop criterion variables
        rk = V1 * D_inv * U1' + V2 * Z * U1' - Ek
        sk = rho * V2' * (Ek - Ekm) * U1
        matrix_norms = [norm(Ek), norm(V2 * Z * U1'), norm(V1 * D_inv * U1')]
        primal_upper_bound = eps_abs * sqrt(m*n) + eps_rel * maximum(matrix_norms)
        aux_var = norm(V2' * Lambda * U1)
        dual_upper_bound = eps_abs * sqrt((n-r)*r) + eps_rel * rho * aux_var
        # Checks stop criterion
        if (norm(rk) <= primal_upper_bound) && (norm(sk) <= dual_upper_bound)
            break
        end
        # Makes Ek the new Ek-1
        Ekm = Ek
        elapsed_time = time() - start_time
        if elapsed_time > time_limit
            println("TimeLimit: ADMM exceed time limit to solve the problem.")
            return "-"
        end
    end
    # Calculates output matrix H
    H = V1 * D_inv * U1' + V2 * Z * U1'
    return H
end

function admm_p1(A::Matrix{Float64})
    m, n = size(A)
    V = 2 * (I(m*n) + kron(A * A', A' * A))
    LLt = cholesky(V)
    H = rand(n, m)
    E = rand(n, m)
    Lambda = rand(n, m)
    Gamma = rand(m, n)
    rho = 4.0
    iter = 1000
    for k in 1:iter
        # Updates H
        Z = 2 * (E - Lambda + A' * (A + Gamma) * A')
        if any(isnan, Z)
            println("Z has nan.")
            break
        end
        z = vec(Z)
        y = LLt.L \ z
        h = LLt.U \ y
        if norm(h - hz) > 10^(-5)
            println("Linear system fail: $k")
            break
        end
        H = reshape(h, size(H))
        if any(isnan, H)
            println("H has nan.")
            break
        end
        # Updates E
        Y = H + Lambda
        E = soft_thresholding_matrix(Y, 1/rho)
        if any(isnan, E)
            println("E has nan.")
            break
        end
        # Updates Lambda
        Lambda = Lambda + H - E
        if any(isnan, Lambda)
            println("Lambda has nan.")
            break
        end
        # Updates Gamma
        Gamma = Gamma + A - A * H * A
        if any(isnan, Gamma)
            println("Gamma has nan.")
            break
        end
        # Print
        # if k % 100 == 0
        #     println("Iteration: $k")
        #     # println("Primal residual: $(norm(rp))")
        #     # println("Dual residual: $(norm(rd))")
        #     println("Obj. Func. Val.: $(norm(H, 1))")
        # end
    end
    return H
end