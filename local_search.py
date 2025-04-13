import matlab.engine
import os

from utility import read_matrix, get_m_n_r_d_idx_from_matrix_filepath

def local_search_procedure(experiment, matrix_filepath, func_name, hatA_flag):
    if func_name not in ["LSFI_Det", "LSFI_Det_Symmetric", "LSFI_Det_P3"]:
        raise Exception("Error: Function does not exist.")

    # Initiates Matlab engine
    eng = matlab.engine.start_matlab()

    # Gets current directory
    current_directory = os.getcwd()

    # Changes Matlab working directory
    eng.cd(current_directory, nargout=0)

    # Reads matrix
    A = read_matrix(matrix_filepath)

    # Gets matrix data
    m, n, r, d, idx = get_m_n_r_d_idx_from_matrix_filepath(matrix_filepath)

    # Specifies LS matrix name
    LS_matrix_name = matrix_filepath.split("/")[-1].replace("matrix", "LS_matrix")

    # Specifies path to save LS matrix
    save_path = f"Local_Search_Matrices/Experiment_{experiment}/{LS_matrix_name}"

    # Calls Local Search Procedure
    eng.call_local_search_procedure(matrix_filepath, r, m, n, func_name, save_path, hatA_flag, nargout=0)

    # Reads LS matrix
    LS_matrix = read_matrix(save_path)

    # Ends Matlab engine
    eng.quit()

    return LS_matrix