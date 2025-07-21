using LinearAlgebra

epsilon = 10^(-5)

function gen_random_rank_r_matrix(m::Int64, n::Int64, r::Int64)
    A = rand(m, n)
    U, S, V = svd(A)
    S_diag = Diagonal(S)
    for i in r+1:min(m, n)
        S_diag[i, i] = 0
    end
    A = U*S_diag*V'
    return A
end

# function matrix_norm_0(A::Matrix{Float64})
#     norm_0 = count(x -> abs(x) > epsilon, A)
#     return norm_0
# end

function matrix_norm_0(A::Union{Matrix{Float64}, Matrix{Int64}})
    norm_0 = 0
    for a in A
        if abs(a) > epsilon
            norm_0 += 1
        end
    end
    return norm_0
end

function calculate_rank(A::Union{Matrix{Float64}, Matrix{Int64}})
    U, S, V = svd(A)
    S = Diagonal(S)
    rank = 0
    for s in S
        if abs(s) > epsilon
            rank += 1
        end
    end
    return rank
end

function count_singular_values(S::Diagonal)
    rank = 0
    for a in S
        if abs(a) > epsilon
            rank += 1
        end
    end
    return rank
end