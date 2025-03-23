using Gurobi, LinearAlgebra
using JuMP

function gurobi_solver(data::DataInst, constraints::Vector{String}, opt_tol::Float64=10^(-8))
    A = data.A
    U,S,V = svd(A);
    Ur = U[:,1:data.r]; Sr = S[1:data.r]; Vr = V[:,1:data.r];
    model = Model(Gurobi.Optimizer)
    @variable(model, t)
    if "Sym" in constraints
        @variable(model, H[1:data.n, 1:data.m],Symmetric)
    else
        @variable(model, H[1:data.n, 1:data.m])
   end
    @objective(model, Min, t)
    @constraint(model, [t; vec(H)] in MOI.NormOneCone(1 + data.n*data.m))
    @constraint(model,-1e-6 .<= Vr'*H*Ur - diagm(1 ./Sr).<= 1e-6)
    set_attribute(model, "OptimalityTol", 1e-5)
    set_attribute(model, "TimeLimit", 180)

    optimize!(model)

    status = termination_status(model)
    if status == MOI.OPTIMAL
        H_star = value.(H)
        return H_star
    else
        throw(ErrorException("Model was not optimized successfully. Status Code: $status"))
    end
end