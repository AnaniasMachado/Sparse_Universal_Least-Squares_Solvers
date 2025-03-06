using MATLAB

function local_search_procedure()
# Initiates a Matlab session
matlab = MATLAB.Matlab()

# Adds directory path to Matlab path
push!(matlab.path, "../")

# Calls Matlab function
resultado = matlab.minha_funcao(arg1, arg2)

# Closes the Matlab session
close(matlab)