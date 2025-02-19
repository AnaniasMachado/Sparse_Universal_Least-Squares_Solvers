function call_local_search_procedure(A_path, r, m, n, R, C, funcName, savePath)
    % Loads matrix A in the specified path
    A = load(A_path);
    
    % Checks if A is a 2D matrix
    if isfield(A, 'matrix')  % Checks if structure contains matrix
        A = A.A;  % Atributes matrix if in structure
    end
    
    % Calls an appropriate function with given data parameters
    switch funcName
        case 'LSFI_Det'
            [normLSFI_Det, timeLSFI_Det, swaps_LSFI_Det, R_LSFI_Det, C_LSFI_Det] = LSFI_Det(A, r, m, n, R, C);
            resultMatrix = inv(A(R_LSFI_Det, C_LSFI_Det));  % Stores resulting matrix
        case 'LSFI_Det_Symmetric'
            [norm, time, swaps, R] = LSFI_Det_Symmetric(A, r, m, R);
            resultMatrix = inv(A(R, R));  % Stores resulting matrix
        case 'LSFI_Det_P3'
            [norm, time, swaps, C_out] = LSFI_Det_P3(A, r, n, R, C);
            A_hat = A(:, C_out);
            resultMatrix = (A_hat' * A_hat) \ A_hat';  % Stores resulting matrix
        otherwise
            error('Function not recognized. Use "LSFI_Det", "LSFI_Det_Symmetric" or "LSFI_Det_P3".');
    end
    
    % Saves resulting matrix in the specified path
    save(savePath, 'matrix');
end