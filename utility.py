import numpy as np

epsilon = 10 ** -5

def generate_random_rank_r_matrix(m):
    n = int(0.5 * m)
    r = int(np.floor(0.25 * m))
    while True:
        A = np.random.rand(m, n)
        U, S, VT = np.linalg.svd(A)
        S_bar = np.zeros((m, n))
        if S.shape[0] < r:
            continue
        for i in range(r):
            S_bar[i, i] = S[i]
        A = np.dot(U, np.dot(S_bar, VT))
        return A

def generate_random_rank_r_square_matrix_vector_b(m):
    r = int(np.floor(0.5 * m))
    while True:
        A = np.random.rand(m, m)
        U, S, VT = np.linalg.svd(A)
        S_bar = np.zeros((m, m))
        if S.shape[0] < r:
            continue
        for i in range(r):
            S_bar[i, i] = S[i]
        A = np.dot(U, np.dot(S_bar, VT))
        coefficients = np.random.rand(A.shape[1])
        b = np.dot(A, coefficients)
        return A, b

def vec_1_norm(x):
    return np.linalg.norm(x, ord=1)

def matrix_vec_1_norm(H):
    return np.linalg.norm(H.flatten(), ord=1)

def matrix_vec_0_norm(H):
    return np.linalg.norm(H.flatten(), ord=0)

def matrix_frobenius_norm(H):
    return np.linalg.norm(H, ord="fro")

def matrix_rank(H):
    U, S, VT = np.linalg.svd(H)
    rank = 0
    for i in range(S.shape[0]):
        if S[i] > epsilon:
            rank += 1
    return rank

def calculate_problem_results(A, H, problem):
    results = dict()

    results[f"{problem}_||H||_1"] = matrix_vec_1_norm(H)
    results[f"{problem}_||H||_0"] = matrix_vec_0_norm(H)
    results[f"{problem}_r(H)"] = matrix_rank(H)
    AHA = np.dot(A, np.dot(H, A))
    results[f"{problem}_||AHA - A||_F"] = matrix_frobenius_norm(AHA - A)
    HAH = np.dot(H, np.dot(A, H))
    results[f"{problem}_||HAH - H||_F"] = matrix_frobenius_norm(HAH - H)
    AH_T = np.dot(A, H).T
    AH = np.dot(A, H)
    results[f"{problem}_||(AH)^T - AH||_F"] = matrix_frobenius_norm(AH_T - AH)
    ATAH = np.dot(A.T, np.dot(A, H))
    results[f"{problem}_||A^TAH - AT||_F"] = matrix_frobenius_norm(ATAH - A.T)
    return results

def calculate_problem_results_3(A, H, b, problem):
    results = dict()

    results[f"{problem}_||H||_1"] = matrix_vec_1_norm(H)
    results[f"{problem}_||H||_0"] = matrix_vec_0_norm(H)
    results[f"{problem}_r(H)"] = matrix_rank(H)
    AHA = np.dot(A, np.dot(H, A))
    results[f"{problem}_||AHA - A||_F"] = matrix_frobenius_norm(AHA - A)
    HAH = np.dot(H, np.dot(A, H))
    results[f"{problem}_||HAH - H||_F"] = matrix_frobenius_norm(HAH - H)
    HA_T = np.dot(H, A).T
    HA = np.dot(H, A)
    results[f"{problem}_||(HA)^T - HA||_F"] = matrix_frobenius_norm(HA_T - HA)
    AHb = np.dot(A, np.dot(H, b))
    results[f"{problem}_||AHb - b||_1"] = vec_1_norm(AHb - b)
    return results

def problem_1_norm_P1_viable_solution(A, H, m, n):
    AHA = np.dot(A, np.dot(H, A))
    for i in range(m):
        for j in range(n):
            if np.abs(AHA[i, j] - A[i, j]) > epsilon:
                return False
    return True

def problem_1_norm_PLS_viable_solution(A, H, m, n):
    ATAH = np.dot(A.T, np.dot(A, H))
    AT = A.T
    for i in range(n):
        for j in range(m):
            if np.abs(ATAH[i, j] - AT[i, j]) > epsilon:
                return False
    return True

def problem_1_norm_P1_PLS_viable_solution(A, H, m, n):
    AHA = np.dot(A, np.dot(H, A))
    for i in range(m):
        for j in range(n):
            if np.abs(AHA[i, j] - A[i, j]) > epsilon:
                return False
    ATAH = np.dot(A.T, np.dot(A, H))
    AT = A.T
    for i in range(n):
        for j in range(m):
            if np.abs(ATAH[i, j] - AT[i, j]) > epsilon:
                return False
    return True

def problem_1_norm_P1_P3_viable_solution(A, H, m, n):
    AHA = np.dot(A, np.dot(H, A))
    for i in range(m):
        for j in range(n):
            if np.abs(AHA[i, j] - A[i, j]) > epsilon:
                return False
    AH_T = np.dot(A, H).T
    AH = np.dot(A, H)
    for i in range(m):
        for j in range(m):
            if np.abs(AH_T[i, j] - AH[i, j]) > epsilon:
                return False
    return True

def problem_1_norm_PN_viable_solution(A, H, b, m):
    epsilon = 10 ** -5
    AHb = np.dot(A, np.dot(H, b))
    for i in range(m):
        if np.abs(AHb[i] - b[i]) > epsilon:
            print("AHb - b == 0 violation")
            return False
    return True

def problem_1_norm_P1_P4_viable_solution(A, H, m):
    epsilon = 10 ** -5
    AHA = np.dot(A, np.dot(H, A))
    for i in range(m):
        for j in range(m):
            if np.abs(AHA[i, j] - A[i, j]) > epsilon:
                return False
    HA_T = np.dot(H, A).T
    HA = np.dot(H, A)
    for i in range(m):
        for j in range(m):
            if np.abs(HA_T[i, j] - HA[i, j]) > epsilon:
                return False
    return True

def problem_1_norm_MSN_viable_solution(A, H, m):
    epsilon = 10 ** -5
    HAA_T = np.dot(H, np.dot(A, A.T))
    A_T = A.T
    for i in range(m):
        for j in range(m):
            if np.abs(HAA_T[i, j] - A_T[i, j]) > epsilon:
                return False
    return True