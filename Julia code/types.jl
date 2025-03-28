using JuMP

struct DataInst
    A::Matrix{Float64}
    m::Int64
    n::Int64
    r::Int64
end

mutable struct GurobiInst
    model::JuMP.Model
    H::Matrix{VariableRef}
end

# struct DRSProjData
#     R::Matrix{Float64}
#     S::Matrix{Float64}
#     T::Matrix{Float64}
#     RMTSM::Matrix{Float64}
#     SSMP::Matrix{Float64}
#     RMPR::Matrix{Float64}
#     VG::Matrix{Float64}
# end

struct DRSProjDataSimple
    R::Matrix{Float64}
    S::Matrix{Float64}
    T::Matrix{Float64}
    RMP::Matrix{Float64}
    SMP::Matrix{Float64}
    T_factor::Matrix{Float64}
    U::Matrix{Float64}
end

struct DRSProjDataP123
    AMP::Matrix{Float64}
end