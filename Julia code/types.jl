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

# mutable struct GurobiInst
#     model::JuMP.Model
#     H::JuMP.Matrix{VariableRef}
# end