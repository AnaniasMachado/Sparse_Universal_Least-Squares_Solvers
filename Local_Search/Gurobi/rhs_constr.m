% Right-hand side constraint construction (rhs_constr.m)

% DESCRIPTION: Given which constraints are included, generate rhs

% INPUT: nr - # of rows in constraint matrix
%        set_prop - cell array, contains status of constraint type (incl)
%        L/U - box constraints for H (lower/upper bounds)

function[z_rhs] = rhs_constr(m,n,A,set_prop,L,M)
epsilon=10^(-4);
z_rhs = [];

% 1-norm objective reformulation: H + Z >= 0 and Z - H >= 0
z_rhs = [z_rhs, zeros(1,2*n*m)];
        

% P1: AHA = A, rhs = A(i,j) for i=1:m, j=1:n
if(strcmp(set_prop{1},'included'))
    
    % SVD formulation: r*r entries [use when rank(A) << min(n,m)]
    
    [U,~,V] = svd(A,'econ');
    r = rank(A);
    
    for i=1:r
        for j=1:r
            
            rhs_SVD = 0; 
            
            for k=1:m
                for l=1:n
                    
                    rhs_SVD = rhs_SVD + U(k,i)*A(k,l)*V(l,j);
                    
                end % end of for loop, for l=1:n
            end % end of for loop, for k=1:m
                                
            z_rhs = [z_rhs, rhs_SVD];
            
        end % end of for loop, for j=1:r
    end % end of for loop, for i=1:r
    
%     % All right hand sides A(i,j) for all i,j (n*m entries)
%     
%     for i=1:m
%         for j=1:n
%             
%             z_rhs = [z_rhs, A(i,j)];
%             
%         end
%     end

end % end of if statement


% P3: (AH) - (AH)' = 0, rhs = 0 
if(strcmp(set_prop{2},'included'))
    
    % All right hand sides are 0 (m*m entries)
    
    z_rhs = [z_rhs, zeros(1,(m*m - m)/2)];
    
end % end of if statement


% P4: (HA) - (HA)' = 0, rhs = 0
if(strcmp(set_prop{3},'included'))
    
    % All right hand sides are 0 (n*n entries)
    
    z_rhs = [z_rhs, zeros(1,(n*n - n)/2)];
    
end % end of if statement


% Nonsym P2: HAH - H = 0, rhs = 0
if(strcmp(set_prop{4},'included'))
    
    % All right hand sides are 0 (n*m entries)
    
    z_rhs = [z_rhs, zeros(1,n*m)];
    
    
end % end of if statement


% Prod. Tightening of P1: <A,K[j,i]> - H(j,i) = 0; H(j,i)(AHA - A) = 0, rhs = 0
if(strcmp(set_prop{5},'included'))
    
    % All right hand sides 0 ([n*m]^2 entries)
    
    z_rhs = [z_rhs, zeros(1,(n*m) + (n*m)^2)];

end % end of if statement    
    
    
% McCormick Lifting Ineq: 
if(strcmp(set_prop{6},'included'))
    
    % Type 1: -L(i,k)*L(l,j) rhs
    
    for i=1:n
        for j=1:m
            
            for k=1:m
                for l=1:n
                    
                    z_rhs = [z_rhs, -L(i,k)*L(l,j)];
                    
                end % end of for loop, for l=1:n
            end % end of for loop, for k=1:m
        end % end of for loop, for j=1:m
    end % end of for loop, for i=1:n
    
    % Type 2: L(i,k)*M(l,j) rhs
    
    for i=1:n
        for j=1:m
            
            for k=1:m
                for l=1:n
                    
                    z_rhs = [z_rhs, L(i,k)*M(l,j)];
                    
                end % end of for loop, for l=1:n
            end % end of for loop, for k=1:m
        end % end of for loop, for j=1:m
    end % end of for loop, for i=1:n
    
    % Type 3: -M(i,k)*M(l,j) rhs
    
    for i=1:n
        for j=1:m
            
            for k=1:m
                for l=1:n
                    
                    z_rhs = [z_rhs, -M(i,k)*M(l,j)];
                    
                end % end of for loop, for l=1:n
            end % end of for loop, for k=1:m
        end % end of for loop, for j=1:m
    end % end of for loop, for i=1:n
    
end % end of if statement


% Quadratic Lifting Ineq:
if(strcmp(set_prop{7},'included'))
    
    % Only for secant inequalities and upper/lower bounds on t_p{ij}
    
    % Secant inequalities:
    
    for p=1:2
        
        for i=1:n
            for j=1:m
                
                ij = (i-1)*m + j;
                
                sec_rhs = (a(p,ij)^(2)) - a(p,ij)*((b(p,ij)^2 - a(p,ij)^2)/(b(p,ij) - a(p,ij)));
                
                z_rhs = [z_rhs, sec_rhs];
                
            end % end of for loop, for j=1:m
        end % end of for loop, for i=1:n
        
    end % end of for loop, for p=1:2
    
    % Upper/Lower bounds on t_p{ij} 
    
    % Upper bound: t_p{ij} <= b(p,ij)
    for p=1:2
        
        for i=1:n
            for j=1:m
                
                ij = (i-1)*m + j;
                
                z_rhs = [z_rhs,b(p,ij)];
                
            end % end of for loop, for j=1:m
        end % end of for loop, for i=1:n
        
    end % end of for loop, for p=1:2
    
    % Lower bound: -t_p{ij} <= -a(p,ij)
    for p=1:2
        
        for i=1:n
            for j=1:m
                
                ij = (i-1)*m + j;
                
                z_rhs = [z_rhs,-a(p,ij)];
                
            end % end of for loop, for j=1:m
        end % end of for loop, for i=1:n
        
    end % end of for loop, for p=1:2
                  
    
end % end of if statement

% Box constraints:
if(strcmp(set_prop{8},'included'))
    
    % Box constraints: Lower bounds (L), then upper bounds (M)
    
    for i=1:n
        for j=1:m
            
            % -H(i,j) <= -L(i,j)
            z_rhs = [z_rhs, -L(i,j)];
            
        end % end of for loop, for j=1:m
    end % end of for loop, for i=1:n
    
    for i=1:n
        for j=1:m
            
            % H(i,j) <= M(i,j)
            z_rhs = [z_rhs, M(i,j)];
            
        end % end of for loop, for j=1:m
    end % end of for loop, for i=1:n
    
end % end of if statement

% Property 10: Symmetric H - H' =0 ,rhs =0

if(strcmp(set_prop{10},'included'))
    
    % All right hand sides are 0 ((m*m -m)/2 entries)
    
    %%%z_rhs = [z_rhs, zeros(1,(m*m - m)/2)];
    z_rhs = [z_rhs, -epsilon*ones(1,(m*m - m))];
    
end 

% Property 11: P2L HAA+ - H =0 ,rhs =0

if(strcmp(set_prop{11},'included'))
    
    % All right hand sides are 0 (n*m entries)
    z_rhs = [z_rhs, zeros(1,n*m)];
    
end 
end
    



