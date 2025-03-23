import numpy as np
import os
import scipy.io
import scipy.sparse
import re

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

# def matrix_vec_0_norm(H):
#     return np.linalg.norm(H.flatten(), ord=0)

def matrix_vec_0_norm(H):
    total = 0
    for x in H.flatten():
        if x > epsilon:
            total += 1
    return total

def matrix_vec_inf_norm(H):
    return max([np.abs(h) for h in H.flatten()])

def matrix_frobenius_norm(H):
    return np.linalg.norm(H, ord="fro")

def matrix_2_1_norm(H):
    total = 0
    for i in range(H.shape[0]):
        total += np.linalg.norm(H[i, :], ord=2)
    return total

def matrix_2_0_norm(H):
    total = 0
    for i in range(H.shape[0]):
        if np.linalg.norm(H[i, :], ord=2) > epsilon:
            total += 1
    return total

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

def calculate_problem_results_5(A, H, problem):
    results = dict()

    results[f"{problem}_r(H)"] = matrix_rank(H)
    results[f"{problem}_||H||_1"] = matrix_vec_1_norm(H)
    results[f"{problem}_||H||_0"] = matrix_vec_0_norm(H)
    results[f"{problem}_||H||_2,1"] = matrix_2_1_norm(H)
    results[f"{problem}_||H||_2,0"] = matrix_2_0_norm(H)
    return results

def calculate_problem_results_6(A, H, problem):
    results = dict()

    results[f"{problem}_||H||_0"] = matrix_vec_0_norm(H)
    results[f"{problem}_||H||_1"] = matrix_vec_1_norm(H)
    return results

def save_H(experiment, problem, m, H):
    filepath = f"./H_stars/experiment_{experiment}_problem_{problem}_m_{m}.npy"
    np.save(filepath, H)

def get_experiment_matrices_filepath(experiment):
    # Defines directory path
    diretory = f"./Experiment_Matrices/Experiment_{experiment}"

    # Creates a list with every matrix name
    matrices_filepath = []
    for filename in os.listdir(diretory):
        filepath = os.path.join(diretory, filename)
        if os.path.isfile(filepath):
            matrix_filepath = f"{diretory}/{filename}"
            matrices_filepath.append(matrix_filepath)
    return matrices_filepath

def read_matrix(matrix_filepath):
    # Loads file .mat
    mat = scipy.io.loadmat(matrix_filepath)

    # Reads matrix
    matrix = mat['matrix']

    if scipy.sparse.issparse(matrix):
        # Converts to dense matrix
        dense_matrix = matrix.toarray()
        return dense_matrix
    return matrix

def get_m_n_r_d_from_matrix_filepath(matrix_filepath):
    # Using regex to find values of m, n, r and d
    match = re.search(r'm(\d+)_n(\d+)_r(\d+)_d(\d+)', matrix_filepath)

    if match:
        m = int(match.group(1))
        n = int(match.group(2))
        r = int(match.group(3))
        d = int(match.group(4))
        
        return [m, n, r, d / 100]
    else:
        raise Exception("ReadError: Could not read all values.")

def problem_1_norm_P1_viable_solution(A, H, m, n):
    m, n = A.shape
    AHA = np.dot(A, np.dot(H, A))
    for i in range(m):
        for j in range(n):
            if np.abs(AHA[i, j] - A[i, j]) > epsilon:
                return False
    return True

def problem_1_norm_PLS_viable_solution(A, H, m, n):
    m, n = A.shape
    ATAH = np.dot(A.T, np.dot(A, H))
    AT = A.T
    for i in range(n):
        for j in range(m):
            if np.abs(ATAH[i, j] - AT[i, j]) > epsilon:
                return False
    return True

def problem_1_norm_P1_PLS_viable_solution(A, H, m, n):
    m, n = A.shape
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
    m, n = A.shape
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
    m, n = A.shape
    AHb = np.dot(A, np.dot(H, b))
    for i in range(m):
        if np.abs(AHb[i] - b[i]) > epsilon:
            print("AHb - b == 0 violation")
            return False
    return True

def problem_1_norm_P1_P4_viable_solution(A, H, m):
    m, n = A.shape
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
    m, n = A.shape
    epsilon = 10 ** -5
    HAA_T = np.dot(H, np.dot(A, A.T))
    A_T = A.T
    for i in range(m):
        for j in range(m):
            if np.abs(HAA_T[i, j] - A_T[i, j]) > epsilon:
                return False
    return True

def problem_1_norm_P1_P3_P4_viable_solution(A, H, m, n):
    m, n = A.shape
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
    HA_T = np.dot(H, A).T
    HA = np.dot(H, A)
    for i in range(n):
        for j in range(n):
            if np.abs(HA_T[i, j] - HA[i, j]) > epsilon:
                return False
    return True

def problem_1_norm_PLS_P4_viable_solution(A, H, m, n):
    m, n = A.shape
    ATAH = np.dot(A.T, np.dot(A, H))
    A_T = A.T
    for i in range(n):
        for j in range(m):
            if np.abs(ATAH[i, j] - A_T[i, j]) > epsilon:
                return False
    HA_T = np.dot(H, A).T
    HA = np.dot(H, A)
    for i in range(n):
        for j in range(n):
            if np.abs(HA_T[i, j] - HA[i, j]) > epsilon:
                return False
    return True

def problem_1_norm_PMN_P3_viable_solution(A, H, m, n):
    m, n = A.shape
    AATHT = np.dot(A, np.dot(A.T, H.T))
    for i in range(m):
        for j in range(n):
            if np.abs(AATHT[i, j] - A[i, j]) > epsilon:
                return False
    AH_T = np.dot(A, H).T
    AH = np.dot(A, H)
    for i in range(m):
        for j in range(m):
            if np.abs(AH_T[i, j] - AH[i, j]) > epsilon:
                return False
    return True

def problem_1_norm_PMX_viable_solution(A, H, m, n):
    m, n = A.shape
    LeftHS = np.dot(A, np.dot(A.T, H.T)) + np.dot(H.T, np.dot(A.T, A))
    RightHS = 2*A
    for i in range(m):
        for j in range(n):
            if np.abs(LeftHS[i, j] - RightHS[i, j]) > epsilon:
                return False
    return True

def problem_1_norm_P1_sym_viable_solution(A, H, m, n):
    m, n = A.shape
    AHA = np.dot(A, np.dot(H, A))
    for i in range(m):
        for j in range(n):
            if np.abs(AHA[i, j] - A[i, j]) > epsilon:
                return False
    H_T = H.T
    for i in range(n):
        for j in range(m):
            if np.abs(H[i, j] - H_T[i, j]) > epsilon:
                return False
    return True

def problem_1_norm_P1_P2_P3_viable_solution(A, H, m, n):
    m, n = A.shape
    AHA = np.dot(A, np.dot(H, A))
    for i in range(m):
        for j in range(n):
            if np.abs(AHA[i, j] - A[i, j]) > epsilon:
                return False
    HAH = np.dot(H, np.dot(A, H))
    for i in range(n):
        for j in range(m):
            if np.abs(HAH[i, j] - H[i, j]) > epsilon:
                return False
    AH_T = np.dot(A, H).T
    AH = np.dot(A, H)
    for i in range(m):
        for j in range(m):
            if np.abs(AH_T[i, j] - AH[i, j]) > epsilon:
                return False
    return True