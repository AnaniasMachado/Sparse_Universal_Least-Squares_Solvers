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

constraints_dict = Dict(
    "P1" => add_constraint_P1,
    "P3" => add_constraint_P3,
    "P4" => add_constraint_P4,
    "PLS" => add_constraint_PLS,
    "PMN" => add_constraint_PMN,
    "PMX" => add_constraint_PMX
)

function add_constraints(inst::GurobiInst, data::DataInst, constraints::Vector{String})
    for constraint in constraints
        if constraint != "Sym"
            constraints_dict[constraint](inst, data)
        end
    end
end

function gurobi_solver(data::DataInst, constraints::Vector{String}, opt_tol::Float64=10^(-8))
    null_matrix = zeros(data.n, data.m)

    model = Model(Gurobi.Optimizer)

    if "Sym" in constraints
        @variable(model, H[1:data.n, 1:data.m], Symmetric, lower_bound=-Inf, upper_bound=Inf)
        @variable(model, Z[1:data.n, 1:data.m], Symmetric, lower_bound=-Inf, upper_bound=Inf)
    else
        @variable(model, H[1:data.n, 1:data.m], lower_bound=-Inf, upper_bound=Inf)
        @variable(model, Z[1:data.n, 1:data.m], lower_bound=-Inf, upper_bound=Inf)
    end

    @objective(model, Min, sum(Z[i, j] for i in 1:data.n, j in 1:data.m))

    inst = GurobiInst(model, H)

    add_constraints(inst, data, constraints)

    @constraint(inst.model, Z - H .>= null_matrix, base_name = "Z_minus_H_")
    @constraint(inst.model, Z + H .>= null_matrix, base_name = "Z_plus_H_")

    # set_optimizer_attributes(inst.model, "OptimalityTol" => opt_tol)

    set_attribute(inst.model, "BarConvTol", 1e-5)
    set_attribute(inst.model, "FeasibilityTol", 1e-5)
    set_attribute(inst.model, "OptimalityTol", 1e-5)

    # set_optimizer_attributes(inst.model, "BarConvTol" => 1e-5)
    # set_optimizer_attributes(inst.model, "FeasibilityTol" => 1e-5)
    # set_optimizer_attributes(inst.model, "OptimalityTol" => 1e-5)

    set_optimizer_attribute(inst.model, "TimeLimit", 7200)
    # set_optimizer_attribute(inst.model, "TimeLimit", 600)

    set_optimizer_attribute(inst.model, "LogToConsole", 0)

    # set_optimizer_attribute(inst.model, "LogFile", "gurobi_log.txt")

    optimize!(inst.model)

    status = termination_status(inst.model)
    if status == MOI.OPTIMAL
        H_star = [value(H[i, j]) for i in 1:data.n, j in 1:data.m]
        return H_star
    else
        throw(ErrorException("Model was not optimized successfully. Status Code: $status"))
    end
end