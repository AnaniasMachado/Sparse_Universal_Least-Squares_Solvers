% box_constr_gen.m 

% DESCRIPTION: Generate box constraints for sparse-pseudoinverse H

% INPUT:

% OUTPUT:

function[Ld,Mu] = box_constr_gen(m,n,data_type,data_mat)

if(strcmp(data_type,'provided'))
    
    % Store data_mat in interim matrix D (for ease of syntax)
    D = data_mat;
    
    % Apply some scaling algorithm to generate lower/upper bound matrices
    
    scl = 5; % Scale the entries of D
    eps = 0.1; % Perturb the scaled entries by some epsilon 
    
    Ld = zeros(n,m); % L := lambda
    Mu = zeros(n,m); % M := mu

    for i=1:n
        for j=1:m
            % Identify if D(i,j) is positive or negative
            x = D(i,j);

            if(D(i,j) > 0)
                Ld(i,j) = (-x/scl) - eps;
                Mu(i,j) = ((6*x)/scl) + eps;
            elseif(D(i,j) < 0)
                Ld(i,j) = ((6*x)/scl) - eps;
                Mu(i,j) = (-x/scl) + eps;
            elseif(D(i,j) == 0)
                Ld(i,j) = -eps;
                Mu(i,j) = eps;
            end % end of if/elseif (D(i,j) > || < 0)
        end % end of for loop, for j=1:m
    end % end of for loop, for i=1:n
    
elseif(strcmp(data_type,'generate'))
    
    % Ignore data_mat, which should be an empty array, []
    
    Ld = zeros(n,m); % L := lambda
    Mu = zeros(n,m); % M := mu
    
    % Uniform box constraints
    for i=1:n
        for j=1:m
            
            Ld(i,j) = -10;
            Mu(i,j) = -Ld(i,j);
            
        end
    end
    
    % Non-uniform box constraints (to be constructed)
            
            
    
    
end % end of if/elseif statements, if(strcmp(data_type,[]))