using Gurobi
using JuMP

function gurobi_solver(data::DataInst, constraints::Vector{String}, opt_tol::Float64=10^(-8))

    model = Model(Gurobi.Optimizer)

    if "Sym" in constraints
        @variable(model, H[1:data.n, 1:data.m],Symmetric)
        @variable(model, Z[1:data.n, 1:data.m],Symmetric)

    else
        @variable(model, H[1:data.n, 1:data.m])
        @variable(model, Z[1:data.n, 1:data.m])
   end
    
    A = data.A
   
    @objective(model, Min, sum(Z))

    @constraint(model, Z - H .>= 0)
    @constraint(model, Z + H .>= 0)
    @constraint(model, A*H*A - A .== 0)
    set_attribute(model, "BarConvTol", 1e-5)
    set_attribute(model, "FeasibilityTol", 1e-5)
    set_attribute(model, "OptimalityTol", 1e-5)
    # set_attribute(model, "TimeLimit", 180)
    set_attribute(model, "TimeLimit", 600)

    # set_optimizer_attribute(model, "LogFile", "gurobi_log_3.txt")

    optimize!(model)

    status = termination_status(model)
    if status == MOI.OPTIMAL
        H_star = value.(H)
        return H_star
    else
        throw(ErrorException("Model was not optimized successfully. Status Code: $status"))
    end
end