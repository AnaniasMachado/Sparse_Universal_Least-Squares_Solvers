using Gurobi
using JuMP
using LinearAlgebra

function gurobi_solver_cal(data::DataInst, problem::String, opt_tol::Float64=10^(-5), time_limit::Int64=7200)
    null_matrix = zeros(data.n, data.m)
    U, S, V = svd(data.A, full=true)
    S = Diagonal(S)
    r = count_singular_values(S)
    D = S[1:r, 1:r]
    D_inv = inv(D)
    U1 = U[:, 1:r]
    U2 = U[:, r+1:end]
    V1 = V[:, 1:r]
    V2 = V[:, r+1:end]

    model = Model(Gurobi.Optimizer)

    if problem == "P134"
        @variable(model, W[1:(data.n-r), 1:(data.m-r)], lower_bound=-Inf, upper_bound=Inf)
        @variable(model, T[1:data.n, 1:data.m], lower_bound=-Inf, upper_bound=Inf)

        @objective(model, Min, sum(T[i, j] for i in 1:data.n, j in 1:data.m))

        G = V1 * D_inv * U1'
        H = G + V2 * W * U2'
        @constraint(model, T - H .>= null_matrix, base_name = "T_minus_H_")
        @constraint(model, T + H .>= null_matrix, base_name = "T_plus_H_")

        set_attribute(model, "BarConvTol", 1e-5)
        set_attribute(model, "FeasibilityTol", 1e-5)
        set_attribute(model, "OptimalityTol", opt_tol)

        # set_attribute(inst.model, "DualReductions", 0)

        set_optimizer_attribute(model, "TimeLimit", time_limit)

        set_optimizer_attribute(model, "LogToConsole", 0)

        # set_optimizer_attribute(inst.model, "LogFile", "gurobi_log.txt")

        optimize!(model)

        status = termination_status(model)
        if status == MOI.OPTIMAL
            W_star = [value(W[i, j]) for i in 1:(data.n-r), j in 1:(data.m-r)]
            H_star = G + V2 * W_star * U2'
            return H_star
        else
            throw(ErrorException("Model was not optimized successfully. Status Code: $status"))
        end
    elseif problem == "P13"
        @variable(model, Z[1:(data.n-r), 1:r], lower_bound=-Inf, upper_bound=Inf)
        @variable(model, W[1:(data.n-r), 1:(data.m-r)], lower_bound=-Inf, upper_bound=Inf)
        @variable(model, T[1:data.n, 1:data.m], lower_bound=-Inf, upper_bound=Inf)

        @objective(model, Min, sum(T[i, j] for i in 1:data.n, j in 1:data.m))

        G = V1 * D_inv * U1'
        H = G + V2 * W * U2' + V2 * Z * U1'
        @constraint(model, T - H .>= null_matrix, base_name = "T_minus_H_")
        @constraint(model, T + H .>= null_matrix, base_name = "T_plus_H_")

        set_attribute(model, "BarConvTol", 1e-5)
        set_attribute(model, "FeasibilityTol", 1e-5)
        set_attribute(model, "OptimalityTol", opt_tol)

        # set_attribute(inst.model, "DualReductions", 0)

        set_optimizer_attribute(model, "TimeLimit", time_limit)

        set_optimizer_attribute(model, "LogToConsole", 0)

        # set_optimizer_attribute(inst.model, "LogFile", "gurobi_log.txt")

        optimize!(model)

        status = termination_status(model)
        if status == MOI.OPTIMAL
            Z_star = [value(Z[i, j]) for i in 1:(data.n-r), j in 1:r]
            W_star = [value(W[i, j]) for i in 1:(data.n-r), j in 1:(data.m-r)]
            H_star = G + V2 * W_star * U2' + V2 * Z_star * U1'
            return H_star
        else
            throw(ErrorException("Model was not optimized successfully. Status Code: $status"))
        end
    elseif problem == "P123"
        @variable(model, Z[1:(data.n-r), 1:r], lower_bound=-Inf, upper_bound=Inf)
        @variable(model, T[1:data.n, 1:data.m], lower_bound=-Inf, upper_bound=Inf)

        @objective(model, Min, sum(T[i, j] for i in 1:data.n, j in 1:data.m))

        G = V1 * D_inv * U1'
        H = G + V2 * Z * U1'
        @constraint(model, T - H .>= null_matrix, base_name = "T_minus_H_")
        @constraint(model, T + H .>= null_matrix, base_name = "T_plus_H_")

        set_attribute(model, "BarConvTol", 1e-5)
        set_attribute(model, "FeasibilityTol", 1e-5)
        set_attribute(model, "OptimalityTol", opt_tol)

        # set_attribute(inst.model, "DualReductions", 0)

        set_optimizer_attribute(model, "TimeLimit", time_limit)

        set_optimizer_attribute(model, "LogToConsole", 0)

        # set_optimizer_attribute(inst.model, "LogFile", "gurobi_log.txt")

        optimize!(model)

        status = termination_status(model)
        if status == MOI.OPTIMAL
            Z_star = [value(Z[i, j]) for i in 1:(data.n-r), j in 1:r]
            H_star = G + V2 * Z_star * U1'
            return H_star
        else
            throw(ErrorException("Model was not optimized successfully. Status Code: $status"))
        end
    end
end