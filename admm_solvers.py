import numpy as np

from utility import *

def initial_variables(V1, U1, D_inv, rho):
    Theta = np.dot(V1, U1.T) / matrix_vec_inf_norm(np.dot(V1, U1.T))
    Lambda = Theta / rho
    E = np.dot(V1, np.dot(D_inv, U1.T)) + Lambda
    return Lambda, E

def soft_thresholding(a, kappa):
    if a > kappa:
        return a - kappa
    elif np.abs(a) <= kappa:
        return 0
    elif a < -kappa:
        return a + kappa

def admm1_123(A, rho, epsilon_abs=10 ** -4, epsilon_rel=10 ** -4):
    # Gets dimensions of matrix A
    m = A.shape[0]
    n = A.shape[1]
    # Calculates full singular value decomposition of A
    U, S, VT = np.linalg.svd(A)
    # Calculates rank of A
    r = matrix_rank(A)
    # Calculates variables U1, V1, V2 and D^{-1}
    U1 = U[:, :r]
    V1 = VT.T[:, :r]
    V2 = VT.T[:, r:]
    D_inv = np.zeros((r, r))
    for i in range(r):
            D_inv[i, i] = 1 / S[i]
    # Calculates initial variables
    Lambda, Ekm = initial_variables(V1=V1, U1=U1, D_inv=D_inv, rho=rho)
    Ek = np.zeros((n, m))
    while True:
        # Updates variable J
        J = Ekm - np.dot(V1, np.dot(D_inv, U1.T)) - Lambda
        # Updates variable Zk
        Z = np.dot(V2.T, np.dot(J, U1))
        # Updates variable Y
        Y = np.dot(V1, np.dot(D_inv, U1.T)) + np.dot(V2, np.dot(Z, U1.T)) + Lambda
        # Updates variable Ek
        for i in range(n):
            for j in range(m):
                Ek[i, j] = soft_thresholding(a=Y[i, j], kappa=1/rho)
        # Updates variable Lambdak
        Lambda = Lambda + np.dot(V1, np.dot(D_inv, U1.T)) + np.dot(V2, np.dot(Z, U1.T)) - Ek
        # Calculates stop criterion variables
        rk = np.dot(V1, np.dot(D_inv, U1.T)) + np.dot(V2, np.dot(Z, U1.T)) - Ek
        sk = rho * np.dot(V2.T, np.dot(Ek - Ekm, U1))
        matrix_norms = [
            matrix_frobenius_norm(Ek),
            matrix_frobenius_norm(np.dot(V2, np.dot(Z, U1.T))),
            matrix_frobenius_norm(np.dot(V1, np.dot(D_inv, U1.T)))
        ]
        primal_upper_bound = epsilon_abs * np.sqrt(m*n) + epsilon_rel * max(matrix_norms)
        aux_var = matrix_frobenius_norm(np.dot(V2.T, np.dot(Lambda, U1)))
        dual_upper_bound = epsilon_abs * np.sqrt((n-r)*r) + epsilon_rel * rho * aux_var
        # Checks stop criterion
        if (matrix_frobenius_norm(rk) <= primal_upper_bound) and (matrix_frobenius_norm(sk) <= dual_upper_bound):
            break
        # Makes Ek the new Ek-1
        Ekm = Ek
    # Calculates output matrix H
    H = np.dot(V1, np.dot(D_inv, U1.T)) + np.dot(V2, np.dot(Z, U1.T))
    return H

def admm1e_123(A, rho=3, epsilon=10 ** -4):
    # Gets dimensions of matrix A
    m = A.shape[0]
    n = A.shape[1]
    # Calculates full singular value decomposition of A
    U, S, VT = np.linalg.svd(A)
    # Calculates rank of A
    r = matrix_rank(A)
    # Calculates variables U1, V1, V2 and D^{-1}
    U1 = U[:, :r]
    V1 = VT.T[:, :r]
    V2 = VT.T[:, r:]
    D_inv = np.zeros((r, r))
    for i in range(r):
            D_inv[i, i] = 1 / S[i]
    # Calculates initial variables
    Lambda, Ekm = initial_variables(V1=V1, U1=U1, D_inv=D_inv, rho=rho)
    Ek = np.zeros((n, m))
    while True:
        # Updates variable J
        J = Ekm - np.dot(V1, np.dot(D_inv, U1.T)) - Lambda
        # Updates variable Zk
        Z = np.dot(V2.T, np.dot(J, U1))
        # Updates variable Y
        Y = np.dot(V1, np.dot(D_inv, U1.T)) + np.dot(V2, np.dot(Z, U1.T)) + Lambda
        # Updates variable Ek
        for i in range(n):
            for j in range(m):
                Ek[i, j] = soft_thresholding(a=Y[i, j], kappa=1/rho)
        # Updates variable Lambdak
        Lambda = Lambda + np.dot(V1, np.dot(D_inv, U1.T)) + np.dot(V2, np.dot(Z, U1.T)) - Ek
        # Calculates stop criterion variables
        rk = np.dot(V1, np.dot(D_inv, U1.T)) + np.dot(V2, np.dot(Z, U1.T)) - Ek
        sk = rho * np.dot(V2.T, np.dot(Ek - Ekm, U1))
        # Checks stop criterion
        if (matrix_frobenius_norm(rk) <= epsilon) and (matrix_frobenius_norm(sk) <= epsilon):
            break
        # Makes Ek the new Ek-1
        Ekm = Ek
    # Calculates output matrix H
    H = np.dot(V1, np.dot(D_inv, U1.T)) + np.dot(V2, np.dot(Z, U1.T))
    return H

def admm1_134(A, rho=3, epsilon_abs=10 ** -4, epsilon_rel=10 ** -4):
    # Gets dimensions of matrix A
    m = A.shape[0]
    n = A.shape[1]
    # Calculates full singular value decomposition of A
    U, S, VT = np.linalg.svd(A)
    # Calculates rank of A
    r = matrix_rank(A)
    # Calculates variables U1, U2, V1, V2 and D^{-1}
    U1 = U[:, :r]
    U2 = U[:, r:]
    V1 = VT.T[:, :r]
    V2 = VT.T[:, r:]
    D_inv = np.zeros((r, r))
    for i in range(r):
            D_inv[i, i] = 1 / S[i]
    # Calculates initial variables
    Lambda, Ekm = initial_variables(V1=V1, U1=U1, D_inv=D_inv, rho=rho)
    Ek = np.zeros((n, m))
    while True:
        # Updates variable J
        J = Ekm - np.dot(V1, np.dot(D_inv, U1.T)) - Lambda
        # Updates variable Wk
        W = np.dot(V2.T, np.dot(J, U2))
        # Updates variable Y
        Y = np.dot(V1, np.dot(D_inv, U1.T)) + np.dot(V2, np.dot(W, U2.T)) + Lambda
        # Updates variable Ek
        for i in range(n):
            for j in range(m):
                Ek[i, j] = soft_thresholding(a=Y[i, j], kappa=1/rho)
        # Updates variable Lambdak
        Lambda = Lambda + np.dot(V1, np.dot(D_inv, U1.T)) + np.dot(V2, np.dot(W, U2.T)) - Ek
        # Calculates stop criterion variables
        rk = np.dot(V1, np.dot(D_inv, U1.T)) + np.dot(V2, np.dot(W, U2.T)) - Ek
        sk = rho * np.dot(V2.T, np.dot(Ek - Ekm, U1))
        matrix_norms = [
            matrix_frobenius_norm(Ek),
            matrix_frobenius_norm(np.dot(V2, np.dot(W, U2.T))),
            matrix_frobenius_norm(np.dot(V1, np.dot(D_inv, U1.T)))
        ]
        primal_upper_bound = epsilon_abs * np.sqrt(m*n) + epsilon_rel * max(matrix_norms)
        aux_var = matrix_frobenius_norm(np.dot(V2.T, np.dot(Lambda, U1)))
        dual_upper_bound = epsilon_abs * np.sqrt((n-r)*r) + epsilon_rel * rho * aux_var
        # Checks stop criterion
        if (matrix_frobenius_norm(rk) <= primal_upper_bound) and (matrix_frobenius_norm(sk) <= dual_upper_bound):
            break
        # Makes Ek the new Ek-1
        Ekm = Ek
    # Calculates output matrix H
    H = np.dot(V1, np.dot(D_inv, U1.T)) + np.dot(V2, np.dot(W, U2.T))
    return H

def admm1e_134(A, rho=3, epsilon=10 ** -4):
    # Gets dimensions of matrix A
    m = A.shape[0]
    n = A.shape[1]
    # Calculates full singular value decomposition of A
    U, S, VT = np.linalg.svd(A)
    # Calculates rank of A
    r = matrix_rank(A)
    # Calculates variables U1, U2, V1, V2 and D^{-1}
    U1 = U[:, :r]
    U2 = U[:, r:]
    V1 = VT.T[:, :r]
    V2 = VT.T[:, r:]
    D_inv = np.zeros((r, r))
    for i in range(r):
            D_inv[i, i] = 1 / S[i]
    # Calculates initial variables
    Lambda, Ekm = initial_variables(V1=V1, U1=U1, D_inv=D_inv, rho=rho)
    Ek = np.zeros((n, m))
    while True:
        # Updates variable J
        J = Ekm - np.dot(V1, np.dot(D_inv, U1.T)) - Lambda
        # Updates variable Wk
        W = np.dot(V2.T, np.dot(J, U2))
        # Updates variable Y
        Y = np.dot(V1, np.dot(D_inv, U1.T)) + np.dot(V2, np.dot(W, U2.T)) + Lambda
        # Updates variable Ek
        for i in range(n):
            for j in range(m):
                Ek[i, j] = soft_thresholding(a=Y[i, j], kappa=1/rho)
        # Updates variable Lambdak
        Lambda = Lambda + np.dot(V1, np.dot(D_inv, U1.T)) + np.dot(V2, np.dot(W, U2.T)) - Ek
        # Calculates stop criterion variables
        rk = np.dot(V1, np.dot(D_inv, U1.T)) + np.dot(V2, np.dot(W, U2.T)) - Ek
        sk = rho * np.dot(V2.T, np.dot(Ek - Ekm, U1))
        # Checks stop criterion
        if (matrix_frobenius_norm(rk) <= epsilon) and (matrix_frobenius_norm(sk) <= epsilon):
            break
        # Makes Ek the new Ek-1
        Ekm = Ek
    # Calculates output matrix H
    H = np.dot(V1, np.dot(D_inv, U1.T)) + np.dot(V2, np.dot(W, U2.T))
    return H