using LinearAlgebra

epsilon = 10^(-5)

function gen_random_rank_r_matrix(m::Int64, n::Int64, r::Int64)
    A = rand(m, n)
    U, S, V = svd(A)
    S_diag = Diagonal(S)
    for i in r+1:n
        S_diag[i, i] = 0
    end
    A = U*S_diag*V'
    return A
end

function matrix_norm_0(A::Matrix{Float64})
    norm_0 = count(x -> abs(x) > epsilon, A)
    return norm_0
end

function calculate_rank(A::Matrix{Float64})
    U, S, V = svd(A)
    S = Diagonal(S)
    rank = count(x -> abs(x) > epsilon, S)
    return rank
end

function get_vmsize()
    ru = Vector{RUsage}(undef, 1)
    ccall(:getrusage, Cint, (Cint, Ptr{Cvoid}), 0, ru)
    return ru[1].ru_maxrss
end