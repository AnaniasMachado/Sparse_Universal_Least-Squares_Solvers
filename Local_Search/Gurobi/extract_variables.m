% extract_variables.m (5/2017)

% DESCRIPTION: Given results from Gurobi solution to LP (result.x), and
%   given our knowledge as to how the matrix variables are stacked within result.x, 
%   we can extract from the vector the matrix forms

function[H,Z,MK] = extract_variables(m,n,x,set_prop)

H = zeros(n,m);
Z = zeros(n,m);

for j=1:n
    for i=1:m
        
        index = (j-1)*m + i;
        
        H(j,i) = x(index);
        
        Z(j,i) = x(n*m + index);
        
    end 
end


% \mathcal{K}_{ji} for j=1:n, i=1:m (each K[i,j] is an m x n matrix)
MK = zeros(n*m*(m),n);


if(strcmp(set_prop{4},'included') || strcmp(set_prop{5},'included') || strcmp(set_prop{6},'included') || strcmp(set_prop{7},'included'))

    for j=1:n
        for i=1:m

            % index for matrix in MK, stacked

            %NOTE: K[ij](lk) := h_{ki}h_{jl}
            
            % Index for K[i,j] (an m x n matrix)
            ji = (j-1)*m + i; % mapping index

            for k=1:n
                for l=1:m
                    
                    lk = (k-1)*m + l;

                    if(ji <= lk)

                        % Use variable result from x

                        t_sum = sum(1:ji-1);

                        
                        % UPDTATE(6/6/17): Incorporate Z (nm), w (2nm)
                        MK((ji - 1)*m + l,k) = x(n*m + 2*n*m + (ji)*(n*m) + lk - t_sum);

                    elseif(ji > lk)

                        % Use K[kl](ij) instead of K[ji](lk)
                        
                        kl = (k-1)*m + l; % Index for which matrix (nm total)
                    
                        ij = (j-1)*m + i; % Index for entry in a (m x n) matrix

                        t_sum = sum(1:kl-1);

                        MK((ji - 1)*m + l,k) = x(n*m + 2*n*m + (kl)*(n*m) + ij - t_sum);

                    end % end of if/elseif statement, if(ji >=< kl)

                end % end of inner loop, for l=1:m
            end % end of outer loop, for k=1:n

        end % end of inner loop, for i=1:m
    end % end of outer loop, for j=1:n
    
end % end of if statment, if(strcmp(set_prop{4-7},'included'))
        
