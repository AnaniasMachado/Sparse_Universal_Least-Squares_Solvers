function call_local_search_procedure(A_path, r, m, n, funcName, savePath)
    % Loads matrix A in the specified path
    A = load(A_path);
    
    % Checks if A is a 2D matrix
    if isfield(A, 'matrix')  % Checks if structure contains matrix
        A = A.matrix;  % Atributes matrix if in structure
        A = full(A); % Converts A to a dense matrix
    end

    addpath('Local_Search/nsub');
    [R,C,time] = nsub(A,r);
    
    % Calls an appropriate function with given data parameters
    switch funcName
        case 'LSFI_Det'
            addpath('Local_Search/P1');
            [normLSFI_Det, timeLSFI_Det, swaps_LSFI_Det, R_LSFI_Det, C_LSFI_Det] = LSFI_Det(A, r, m, n, R, C);
            A_inv = inv(A(R_LSFI_Det, C_LSFI_Det));
            matrix = zeros(n, m);
            matrix(C_LSFI_Det, R_LSFI_Det) = A_inv;
        case 'LSFI_Det_Symmetric'
            addpath('Local_Search/P1sym');
            [norm, time, swaps, R] = LSFI_Det_Symmetric(A, r, m, R);
            A_hat = inv(A(R, R));
            matrix = zeros(n, m);
            matrix(R, R) = inv(A_hat);
        case 'LSFI_Det_P3'
            addpath('Local_Search/P13');
            [norm, time, swaps, C_out] = LSFI_Det_P3(A, r, n, R, C);
            A_hat = A(:, C_out);
            matrix_hat = inv(A_hat' * A_hat) * A_hat';
            matrix = zeros(n, m);
            matrix(C_out, :) = matrix_hat;
        otherwise
            error('Function not recognized. Use "LSFI_Det", "LSFI_Det_Symmetric" or "LSFI_Det_P3".');
    end
    
    % Saves resulting matrix in the specified path
    save(savePath, 'matrix');
end