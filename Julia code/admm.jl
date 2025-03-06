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

function admm_p123(A::Matrix{Float64}, rho::Float64, eps_abs::Float64, eps_rel::Float64, fixed_tol::Bool, eps_opt::Float64, time_limit::Int64)
    m, n = size(A)
    U, S, V = svd(A)
    S = Diagonal(S)
    r = count_singular_values(S)
    D = S[1:r, 1:r]
    # D_inv = spdiam(1 ./ D)
    D_inv = inv(D)
    U1 = U[:, 1:r]
    V1 = V[:, 1:r]
    V2 = V[:, r+1:end]
    # D_inv = zeros(r, r)
    # for i in 1:r
    #     D_inv[i, i] = 1 / S[i, i]
    # end

    V1DinvU1T = V1 * D_inv * U1'
    V1DinvU1T_F = norm(V1DinvU1T)

    V2V2T = V2 * V2'
    U1U1T = U1 * U1'

    Z = zeros(n - r, r)
    H = zeros(n, m)
    Lambda, Ekm = variables_initialization(V1, U1, D_inv, rho)
    start_time = time()
    k = 0
    while true
        k += 1
        # Updates J
        # J = Ekm - V1DinvU1T - Lambda
        # Updates Zk
        # Z = V2' * J * U1
        # V2ZU1T = V2 * Z * U1'
        # V2ZU1T = V2 * V2' * (Ekm - V1DinvU1T - Lambda) * U1 * U1'
        # V2ZU1T = Ekm - V1DinvU1T - Lambda
        # V2ZU1T = V2V2T * (Ekm - V1DinvU1T - Lambda) * U1U1T
        V2ZU1T = V2V2T * (Ekm - Lambda) * U1U1T
        H = V1DinvU1T + V2ZU1T
        # Updates Y
        # Y = H + Lambda
        # Updates Ek
        # Ek = soft_thresholding_matrix(Y, 1/rho)
        Ek = soft_thresholding_matrix(H + Lambda, 1/rho)
        res_infeas = H - Ek
        # Updates Lambdak
        # Lambda = Lambda + V1DinvU1T + V2ZU1T - Ek
        Lambda += res_infeas
        # Calculates stop criterion variables
        # rk = V1DinvU1T + V2ZU1T - Ek
        # sk = rho * V2' * (Ek - Ekm) * U1
        rk_F = norm(res_infeas)
        sk_F = rho * norm(V2' * (Ek - Ekm) * U1)
        if k % 100 == 0
            println("ADMM iteration: $k")
            println("ADMM primal residual: $rk_F)")
            println("ADMM dual residual: $sk_F)")
        end
        if !fixed_tol
            matrix_norms = [norm(Ek), norm(V2ZU1T), V1DinvU1T_F]
            primal_upper_bound = eps_abs * sqrt(m*n) + eps_rel * maximum(matrix_norms)
            # aux_var = norm(V2' * Lambda * U1)
            # dual_upper_bound = eps_abs * sqrt((n-r)*r) + eps_rel * rho * aux_var
            dual_upper_bound = eps_abs * sqrt((n-r)*r) + eps_rel * rho * norm(V2' * Lambda * U1)
            # Checks stop criterion
            # if (norm(rk) <= primal_upper_bound) && (norm(sk) <= dual_upper_bound)
            #     H = V1DinvU1T + V2ZU1T
            #     break
            # end
            if (rk_F <= primal_upper_bound) && (sk_F <= dual_upper_bound)
                break
            end
        else
            if (rk_F <= eps_opt) && (sk_F <= eps_opt)
                break
            end
        end
        # Makes Ek the new Ek-1
        Ekm = Ek
        elapsed_time = time() - start_time
        if elapsed_time > time_limit
            println("TimeLimit: ADMM exceed time limit to solve the problem.")
            return "-"
        end
    end
    return H
end

# Gabriel version
# function admm_p123(A::Matrix{Float64}, rho::Float64, eps_abs::Float64, eps_rel::Float64, fixed_tol::Bool, time_limit::Int64)
#     m, n = size(A)
#     U, S, V = svd(A)
#     S = Diagonal(S)
#     r = count_singular_values(S)
#     D = S[1:r, 1:r]
#     D_inv = inv(D)
#     U1 = U[:, 1:r]
#     V1 = V[:, 1:r]
#     V2 = V[:, r+1:end]

#     V1DinvU1T = V1 * D_inv * U1'
#     V1DinvU1T_F = norm(V1DinvU1T)

#     V2V2T = V2 * V2'
#     U1U1T = U1 * U1'

#     Z = zeros(n - r, r)
#     H = zeros(n, m)
#     Lambda, Ekm = variables_initialization(V1, U1, D_inv, rho)
#     Ek = Ekm
#     rho_inv = 1/rho
#     iter=0
#     max_iter = 10^5
#     start_time = time()
#     while true
#         # update Z where W := V2*Z*U1' with J = (- G + E - Λ)
#         W = V2V2T*(Ek - Lambda)*U1U1T
#         # update E
#         H = V1DinvU1T + W
#         # update of E
#         Ek = soft_thresholding_matrix(H + Lambda,rho_inv)
#         # Compute residuals
#         res_infeas = H - Ek
#         primal_res = norm(res_infeas)
#         dual_res = rho*norm(V2'*(Ek-Ekm)*U1)
#         # Update P
#         Lambda += res_infeas
#         # Stopping criteria
#         if !fixed_tol
#             eps_p = sqrt(n*m)*eps_abs + eps_rel*maximum([norm(Ek),norm(W),V1DinvU1T_F])
#             eps_d = sqrt((n-r)*r)*eps_abs + eps_rel*rho*norm(V2'*Lambda*U1)
#         elseif fixed_tol
#             dual_res = rho*norm(V2'*Lambda*U1)
#             eps_p = eps_opt;eps_d = eps_opt
#         else
#             error("stop limit not defined correctly")
#         end
#         iter += 1
#         elapsed_time = time() - start_time
#         if iter == max_iter || elapsed_time >= 7200
#             print("Termination: Max Iter.")
#             break
#         end
#         if (primal_res <= eps_p && dual_res <= eps_d)
#             opt_res = abs(getnorm1(H)-rho*tr(V1DinvU1T'*Lambda))
#             println("Termination: Low Res.")
#             break
#         end
#         Ekm = Ek
#     end
#     return H
# end

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