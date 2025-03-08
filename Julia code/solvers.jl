using Gurobi
using JuMP

# include("types.jl")

function add_constraint_P1(inst::GurobiInst, data::DataInst)
    @constraint(inst.model, data.A * inst.H * data.A .== data.A, base_name = "P1_")
end

function add_constraint_P3(inst::GurobiInst, data::DataInst)
    @constraint(inst.model, data.A * inst.H .== inst.H' * data.A', base_name = "P3_")
end

function add_constraint_P4(inst::GurobiInst, data::DataInst)
    @constraint(inst.model, inst.H * data.A .== data.A' * inst.H', base_name = "P4_")
end

function add_constraint_PLS(inst::GurobiInst, data::DataInst)
    @constraint(inst.model, data.A' * data.A * inst.H .== data.A', base_name = "PLS_")
end

function add_constraint_PMN(inst::GurobiInst, data::DataInst)
    @constraint(inst.model, data.A * data.A' * inst.H' .== data.A, base_name = "PMN_")
end

function add_constraint_PMX(inst::GurobiInst, data::DataInst)
    @constraint(
        inst.model,
        data.A * data.A' * inst.H' + inst.H' * data.A' * data.A .== 2*data.A,
        base_name = "PMX_"
    )
end

function add_constraint_Sym(inst::GurobiInst, data::DataInst)
    @constraint(inst.model, inst.H .== inst.H', base_name = "Sym_")
end

constraints_dict = Dict(
    "P1" => add_constraint_P1,
    "P3" => add_constraint_P3,
    "P4" => add_constraint_P4,
    "PLS" => add_constraint_PLS,
    "PMN" => add_constraint_PMN,
    "PMX" => add_constraint_PMX,
    "Sym" => add_constraint_Sym
)

function add_constraints(inst::GurobiInst, data::DataInst, constraints::Vector{String})
    for constraint in constraints
        constraints_dict[constraint](inst, data)
    end
end

function gurobi_solver(data::DataInst, constraints::Vector{String}, opt_tol::Float64=10^(-8))
    null_matrix = zeros(data.n, data.m)

    model = Model(Gurobi.Optimizer)

    @variable(model, H[1:data.n, 1:data.m])
    @variable(model, Z[1:data.n, 1:data.m])

    @objective(model, Min, sum(Z[i, j] for i in 1:data.n, j in 1:data.m))

    inst = GurobiInst(model, H)

    add_constraints(inst, data, constraints)

    @constraint(inst.model, Z - H .>= null_matrix, base_name = "Z_minus_H_")
    @constraint(inst.model, Z + H .>= null_matrix, base_name = "Z_plus_H_")

    set_optimizer_attributes(model, "OptimalityTol" => opt_tol)

    set_optimizer_attribute(model, "LogToConsole", 0)

    optimize!(model)

    status = termination_status(model)
    if status == MOI.OPTIMAL
        H_star = [value(H[i, j]) for i in 1:data.n, j in 1:data.m]
        return H_star
    else
        throw(ErrorException("Model was not optimized successfully. Status Code: $status"))
    end
end