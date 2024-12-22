import numpy as np
import gurobipy as gp
from gurobipy import GRB

def problem_1_norm_P1_solver(A):
    # Extracts matrix dimensions
    m = A.shape[0]
    n = A.shape[1]

    # Defines auxiliary variables
    J = np.ones((n, m))

    # Creates a new model
    model = gp.Model("lp")

    # Creates variables
    H_var = model.addVars(n, m, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="H")
    Z_var = model.addVars(n, m, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="Z")
    H = []
    Z = []
    for i in range(n):
        H_row = []
        Z_row = []
        for j in range(m):
            H_row.append(H_var[i, j])
            Z_row.append(Z_var[i, j])
        H.append(H_row)
        Z.append(Z_row)
    H = np.array(H)
    Z = np.array(Z)

    # Defines objective function
    objective = np.trace(np.dot(J.T, Z))

    # Sets objective function
    model.setObjective(objective, GRB.MINIMIZE)

    # Adds constraint: P1 given by AHA = A
    AHA = np.dot(A, np.dot(H, A))
    for i in range(m):
        for j in range(n):
            model.addConstr(AHA[i, j] == A[i, j], name=f"constraint_P1_i={i}_j={j}")

    # Adds constraint: Z - H >= 0
    Z_minus_H = Z - H
    for i in range(n):
        for j in range(m):
            model.addConstr(Z_minus_H[i, j] >= 0, name=f"constraint_ineq_1_i={i}_j={j}")

    # Adds constraint: Z + H >= 0
    Z_plus_H = Z + H
    for i in range(n):
        for j in range(m):
            model.addConstr(Z_plus_H[i, j] >= 0, name=f"constraint_ineq_2_i={i}_j={j}")

    # Configures log to not show messages on terminal
    model.setParam("OutputFlag", 0)

    # Defines SoftMemLimit (in GigaBytes)
    # model.setParam('SoftMemLimit', 8)

    # Defines TimeLimit as 2 hours (in seconds)
    model.setParam('TimeLimit', 2*60*60)

    # Optimizes model
    model.optimize()

    # Checks if model was optimized successfully
    if model.status == GRB.OPTIMAL:
        # Extracts model variables
        H_values = np.zeros((n, m))
        for i in range(n):
            for j in range(m):
                H_values[i, j] = H[i, j].X
        return H_values
    else:
        raise RuntimeError(f"Model was not optimized successfully. Status: {model.status}")

def problem_1_norm_PLS_solver(A):
    # Extracts matrix dimensions
    m = A.shape[0]
    n = A.shape[1]

    # Defines auxiliary variables
    J = np.ones((n, m))

    # Creates a new model
    model = gp.Model("lp")

    # Creates variables
    H_var = model.addVars(n, m, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="H")
    Z_var = model.addVars(n, m, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="Z")
    H = []
    Z = []
    for i in range(n):
        H_row = []
        Z_row = []
        for j in range(m):
            H_row.append(H_var[i, j])
            Z_row.append(Z_var[i, j])
        H.append(H_row)
        Z.append(Z_row)
    H = np.array(H)
    Z = np.array(Z)

    # Defines objective function
    objective = np.trace(np.dot(J.T, Z))

    # Sets objective function
    model.setObjective(objective, GRB.MINIMIZE)

    # Adds constraint: PLS given by A^TAH = A^T
    ATAH = np.dot(A.T, np.dot(A, H))
    AT = A.T
    for i in range(n):
        for j in range(m):
            model.addConstr(ATAH[i, j] == AT[i, j], name=f"constraint_PLS_i={i}_j={j}")

    # Adds constraint: Z - H >= 0
    Z_minus_H = Z - H
    for i in range(n):
        for j in range(m):
            model.addConstr(Z_minus_H[i, j] >= 0, name=f"constraint_ineq_1_i={i}_j={j}")

    # Adds constraint: Z + H >= 0
    Z_plus_H = Z + H
    for i in range(n):
        for j in range(m):
            model.addConstr(Z_plus_H[i, j] >= 0, name=f"constraint_ineq_2_i={i}_j={j}")

    # Configures log to not show messages on terminal
    model.setParam("OutputFlag", 0)

    # Defines SoftMemLimit (in GigaBytes)
    # model.setParam('SoftMemLimit', 8)

    # Defines TimeLimit as 2 hours (in seconds)
    model.setParam('TimeLimit', 2*60*60)

    # Optimizes model
    model.optimize()

    # Checks if model was optimized successfully
    if model.status == GRB.OPTIMAL:
        # Extracts model variables
        H_values = np.zeros((n, m))
        for i in range(n):
            for j in range(m):
                H_values[i, j] = H[i, j].X
        return H_values
    else:
        raise RuntimeError(f"Model was not optimized successfully. Status: {model.status}")

def problem_1_norm_P1_PLS_solver(A):
    # Extracts matrix dimensions
    m = A.shape[0]
    n = A.shape[1]

    # Defines auxiliary variables
    J = np.ones((n, m))

    # Creates a new model
    model = gp.Model("lp")

    # Creates variables
    H_var = model.addVars(n, m, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="H")
    Z_var = model.addVars(n, m, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="Z")
    H = []
    Z = []
    for i in range(n):
        H_row = []
        Z_row = []
        for j in range(m):
            H_row.append(H_var[i, j])
            Z_row.append(Z_var[i, j])
        H.append(H_row)
        Z.append(Z_row)
    H = np.array(H)
    Z = np.array(Z)

    # Defines objective function
    objective = np.trace(np.dot(J.T, Z))

    # Sets objective function
    model.setObjective(objective, GRB.MINIMIZE)

    # Adds constraint: P1 given by AHA = A
    AHA = np.dot(A, np.dot(H, A))
    for i in range(m):
        for j in range(n):
            model.addConstr(AHA[i, j] == A[i, j], name=f"constraint_P1_i={i}_j={j}")

    # Adds constraint: PLS given by A^TAH = A^T
    ATAH = np.dot(A.T, np.dot(A, H))
    AT = A.T
    for i in range(n):
        for j in range(m):
            model.addConstr(ATAH[i, j] == AT[i, j], name=f"constraint_PLS_i={i}_j={j}")

    # Adds constraint: Z - H >= 0
    Z_minus_H = Z - H
    for i in range(n):
        for j in range(m):
            model.addConstr(Z_minus_H[i, j] >= 0, name=f"constraint_ineq_1_i={i}_j={j}")

    # Adds constraint: Z + H >= 0
    Z_plus_H = Z + H
    for i in range(n):
        for j in range(m):
            model.addConstr(Z_plus_H[i, j] >= 0, name=f"constraint_ineq_2_i={i}_j={j}")

    # Configures log to not show messages on terminal
    model.setParam("OutputFlag", 0)

    # Defines SoftMemLimit (in GigaBytes)
    # model.setParam('SoftMemLimit', 8)

    # Defines TimeLimit as 2 hours (in seconds)
    model.setParam('TimeLimit', 2*60*60)

    # Optimizes model
    model.optimize()

    # Checks if model was optimized successfully
    if model.status == GRB.OPTIMAL:
        # Extracts model variables
        H_values = np.zeros((n, m))
        for i in range(n):
            for j in range(m):
                H_values[i, j] = H[i, j].X
        return H_values
    else:
        raise RuntimeError(f"Model was not optimized successfully. Status: {model.status}")

def problem_1_norm_P1_P3_solver(A):
    # Extracts matrix dimensions
    m = A.shape[0]
    n = A.shape[1]

    # Defines auxiliary variables
    J = np.ones((n, m))

    # Creates a new model
    model = gp.Model("lp")

    # Creates variables
    H_var = model.addVars(n, m, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="H")
    Z_var = model.addVars(n, m, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="Z")
    H = []
    Z = []
    for i in range(n):
        H_row = []
        Z_row = []
        for j in range(m):
            H_row.append(H_var[i, j])
            Z_row.append(Z_var[i, j])
        H.append(H_row)
        Z.append(Z_row)
    H = np.array(H)
    Z = np.array(Z)

    # Defines objective function
    objective = np.trace(np.dot(J.T, Z))

    # Sets objective function
    model.setObjective(objective, GRB.MINIMIZE)

    # Adds constraint: P1 given by AHA = A
    AHA = np.dot(A, np.dot(H, A))
    for i in range(m):
        for j in range(n):
            model.addConstr(AHA[i, j] == A[i, j], name=f"constraint_P1_i={i}_j={j}")

    # Adds constraint: P3 given by (AH)^T = AH
    AH_T = np.dot(A, H).T
    AH = np.dot(A, H)
    for i in range(m):
        for j in range(m):
            model.addConstr(AH_T[i, j] == AH[i, j], name=f"constraint_P3_i={i}_j={j}")

    # Adds constraint: Z - H >= 0
    Z_minus_H = Z - H
    for i in range(n):
        for j in range(m):
            model.addConstr(Z_minus_H[i, j] >= 0, name=f"constraint_ineq_1_i={i}_j={j}")

    # Adds constraint: Z + H >= 0
    Z_plus_H = Z + H
    for i in range(n):
        for j in range(m):
            model.addConstr(Z_plus_H[i, j] >= 0, name=f"constraint_ineq_2_i={i}_j={j}")

    # Configures log to not show messages on terminal
    model.setParam("OutputFlag", 0)

    # Defines SoftMemLimit (in GigaBytes)
    # model.setParam('SoftMemLimit', 8)

    # Defines TimeLimit as 2 hours (in seconds)
    model.setParam('TimeLimit', 2*60*60)

    # Optimizes model
    model.optimize()

    # Checks if model was optimized successfully
    if model.status == GRB.OPTIMAL:
        # Extracts model variables
        H_values = np.zeros((n, m))
        for i in range(n):
            for j in range(m):
                H_values[i, j] = H[i, j].X
        return H_values
    else:
        raise RuntimeError(f"Model was not optimized successfully. Status: {model.status}")

def problem_2_norm_KKT_solver(A, b):
    # Extracts matrix dimensions
    m = A.shape[0]

    # Creates a new model
    model = gp.Model("nlp")

    # Creates variables
    H_var = model.addVars(m, m, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="H")
    H = []
    for i in range(m):
        H_row = []
        for j in range(m):
            H_row.append(H_var[i, j])
        H.append(H_row)
    H = np.array(H)

    # Defines objective function
    objective = np.dot(np.dot(H, b).T, np.dot(H, b))

    # Sets objective function
    model.setObjective(objective, GRB.MINIMIZE)

    # Adds constraint: AHb = b
    AHb = np.dot(A, np.dot(H, b))
    for i in range(m):
        model.addConstr(AHb[i] == b[i], name=f"constraint_eq_i={i}")

    # Configures log to not show messages on terminal
    model.setParam("OutputFlag", 0)

    # Defines SoftMemLimit (in GigaBytes)
    # model.setParam('SoftMemLimit', 8)

    # Optimizes model
    model.optimize()

    # Checks if model was optimized successfully
    if model.status == GRB.OPTIMAL:
        # Extracts model variables
        lambda_values = []
        for i in range(m):
            lambda_value = model.getConstrByName(f"constraint_eq_i={i}").Pi
            lambda_values.append(lambda_value)
        lambda_values = np.array(lambda_values)
        return lambda_values
    else:
        raise RuntimeError(f"Model was not optimized successfully. Status: {model.status}")

def problem_1_norm_PN_solver(A, b):
    # Solves the following problem to get the lagrange multipliers for the PN equations
    # lambda_values = problem_2_norm_KKT_solver(A=A, b=b)

    # Extracts matrix dimensions
    m = A.shape[0]

    # Defines auxiliary variables
    J = np.ones((m, m))

    # Creates a new model
    model = gp.Model("lp")

    # Creates variables
    H_var = model.addVars(m, m, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="H")
    Z_var = model.addVars(m, m, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="Z")
    H = []
    Z = []
    for i in range(m):
        H_row = []
        Z_row = []
        for j in range(m):
            H_row.append(H_var[i, j])
            Z_row.append(Z_var[i, j])
        H.append(H_row)
        Z.append(Z_row)
    H = np.array(H)
    Z = np.array(Z)

    # Defines objective function
    objective = np.trace(np.dot(J.T, Z))

    # Sets objective function
    model.setObjective(objective, GRB.MINIMIZE)

    # Adds constraint: Primal feasibility AHb = b
    AHb = np.dot(A, np.dot(H, b))
    for i in range(m):
        model.addConstr(AHb[i] == b[i], name=f"constraint_eq_1_i={i}")

    # Adds constraint: Z - H >= 0
    Z_minus_H = Z - H
    for i in range(m):
        for j in range(m):
            model.addConstr(Z_minus_H[i, j] >= 0, name=f"constraint_ineq_1_i={i}_j={j}")

    # Adds constraint: Z + H >= 0
    Z_plus_H = Z + H
    for i in range(m):
        for j in range(m):
            model.addConstr(Z_plus_H[i, j] >= 0, name=f"constraint_ineq_2_i={i}_j={j}")

    # Configures log to not show messages on terminal
    model.setParam("OutputFlag", 0)

    # Defines SoftMemLimit (in GigaBytes)
    # model.setParam('SoftMemLimit', 8)

    # Optimizes model
    model.optimize()

    # Checks if model was optimized successfully
    if model.status == GRB.OPTIMAL:
        # Extracts model variables
        H_values = np.zeros((m, m))
        for i in range(m):
            for j in range(m):
                H_values[i, j] = H[i, j].X
        return H_values
    else:
        raise RuntimeError(f"Model was not optimized successfully. Status: {model.status}")

def problem_1_norm_P1_P4_solver(A, b):
    # Extracts matrix dimensions
    m = A.shape[0]

    # Defines auxiliary variables
    J = np.ones((m, m))

    # Creates a new model
    model = gp.Model("lp")

    # Creates variables
    H_var = model.addVars(m, m, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="H")
    Z_var = model.addVars(m, m, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="Z")
    H = []
    Z = []
    for i in range(m):
        H_row = []
        Z_row = []
        for j in range(m):
            H_row.append(H_var[i, j])
            Z_row.append(Z_var[i, j])
        H.append(H_row)
        Z.append(Z_row)
    H = np.array(H)
    Z = np.array(Z)

    # Defines objective function
    objective = np.trace(np.dot(J.T, Z))

    # Sets objective function
    model.setObjective(objective, GRB.MINIMIZE)

    # Adds constraint: P1 given by AHA = A
    AHA = np.dot(A, np.dot(H, A))
    for i in range(m):
        for j in range(m):
            model.addConstr(AHA[i, j] == A[i, j], name=f"constraint_P1_i={i}_j={j}")

    # Adds constraint: P4 given by (HA)^T = HA
    HA_T = np.dot(H, A).T
    HA = np.dot(H, A)
    for i in range(m):
        for j in range(m):
            model.addConstr(HA_T[i, j] == HA[i, j], name=f"constraint_P3_i={i}_j={j}")

    # Adds constraint: Z - H >= 0
    Z_minus_H = Z - H
    for i in range(m):
        for j in range(m):
            model.addConstr(Z_minus_H[i, j] >= 0, name=f"constraint_ineq_1_i={i}_j={j}")

    # Adds constraint: Z + H >= 0
    Z_plus_H = Z + H
    for i in range(m):
        for j in range(m):
            model.addConstr(Z_plus_H[i, j] >= 0, name=f"constraint_ineq_2_i={i}_j={j}")

    # Configures log to not show messages on terminal
    model.setParam("OutputFlag", 0)

    # Defines SoftMemLimit (in GigaBytes)
    # model.setParam('SoftMemLimit', 8)

    # Optimizes model
    model.optimize()

    # Checks if model was optimized successfully
    if model.status == GRB.OPTIMAL:
        # Extracts model variables
        H_values = np.zeros((m, m))
        for i in range(m):
            for j in range(m):
                H_values[i, j] = H[i, j].X
        return H_values
    else:
        raise RuntimeError(f"Model was not optimized successfully. Status: {model.status}")

def problem_1_norm_MSN_solver(A, b):
    # Extracts matrix dimensions
    m = A.shape[0]

    # Defines auxiliary variables
    J = np.ones((m, m))

    # Creates a new model
    model = gp.Model("lp")

    # Creates variables
    H_var = model.addVars(m, m, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="H")
    Z_var = model.addVars(m, m, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="Z")
    H = []
    Z = []
    for i in range(m):
        H_row = []
        Z_row = []
        for j in range(m):
            H_row.append(H_var[i, j])
            Z_row.append(Z_var[i, j])
        H.append(H_row)
        Z.append(Z_row)
    H = np.array(H)
    Z = np.array(Z)

    # Defines objective function
    objective = np.trace(np.dot(J.T, Z))

    # Sets objective function
    model.setObjective(objective, GRB.MINIMIZE)

    # Adds constraint: MSN given by HAA^T = A^T
    HAA_T = np.dot(H, np.dot(A, A.T))
    A_T = A.T
    for i in range(m):
        for j in range(m):
            model.addConstr(HAA_T[i, j] == A_T[i, j], name=f"constraint_MSN_i={i}_j={j}")

    # Adds constraint: Z - H >= 0
    Z_minus_H = Z - H
    for i in range(m):
        for j in range(m):
            model.addConstr(Z_minus_H[i, j] >= 0, name=f"constraint_ineq_1_i={i}_j={j}")

    # Adds constraint: Z + H >= 0
    Z_plus_H = Z + H
    for i in range(m):
        for j in range(m):
            model.addConstr(Z_plus_H[i, j] >= 0, name=f"constraint_ineq_2_i={i}_j={j}")

    # Configures log to not show messages on terminal
    model.setParam("OutputFlag", 0)

    # Defines SoftMemLimit (in GigaBytes)
    # model.setParam('SoftMemLimit', 8)

    # Optimizes model
    model.optimize()

    # Checks if model was optimized successfully
    if model.status == GRB.OPTIMAL:
        # Extracts model variables
        H_values = np.zeros((m, m))
        for i in range(m):
            for j in range(m):
                H_values[i, j] = H[i, j].X
        return H_values
    else:
        raise RuntimeError(f"Model was not optimized successfully. Status: {model.status}")